---- MODULE StagedGeneCache ----
(* The gene-annotation cache entry: concurrent get-or-fetch callers and a
   clear_cache deleter, as separate OS processes sharing one cache directory.
   StagedCache.tla models the liftover chain over the same staged writer.

   Source modelled (src/pylocuszoom/):
     reference_genes.py::get_genes_for_build  load, fetch on a miss, save
     _gene_cache.py::load_annotations
     _gene_cache.py::save_annotations
     _gene_cache.py::clear_cache
     _http.py::staged_path

   One action is one filesystem call.

     pc        code                                             action
     load      load_annotations: ZipFile(entry)                 GLoad
     fetch     get_genes_for_build: source.fetch                GFetch
     mk        save_annotations: mkdir, staged_path mkstemp     GMk
     open      save_annotations: ZipFile(partial, "w")          GOpen
     write     save_annotations: to_csv, the archive close      GWrite
     replace   staged_path: os.replace                          GReplace
     cleanup   staged_path: unlink(missing_ok)                  GCleanup
     glob      clear_cache: both globs, eagerly                 CGlob
     unlink    clear_cache: cache_file.unlink()                 CUnlink

   ZipFile opens the entry once and reads both members from that descriptor,
   so GLoad observes one content. Anything but a complete archive is a miss
   (the except arms of load_annotations). The two globs are unpacked into one tuple before the loop
   starts, so CGlob is one snapshot and the unlinks follow one at a time. The
   sibling name is ".{entry.name}.XXXX.part"; it ends in ".part", so neither
   "*.csv" nor "annotations_*.zip" matches it and Matches holds the entry only.

   Failure injection: a process in FailProcs raises at FailStage. "fetch"
   propagates to the caller; "create", "write" and "replace" are OSError,
   which save_annotations logs and swallows.

   Ownership: part[n] is the sibling named n with the process that created
   it; holds[p] says p's live staged_path frame owns a sibling.

   No mutex or condition variable exists here, so there is no lock variable,
   wait set or spurious wakeup. Liveness checks that every call returns under
   weak fairness on each process's steps. *)
EXTENDS Naturals, FiniteSets

CONSTANTS Getters, Clearers, InitDest, FailProcs, FailStage

ASSUME /\ Getters \cap Clearers = {}
       /\ "none" \notin Getters /\ "dest" \notin Getters
       /\ InitDest \in {"absent", "good", "corrupt"}
       /\ FailProcs \subseteq Getters
       /\ FailStage \in {"none", "fetch", "create", "write", "replace"}

VARIABLES pc, dest, part, holds, failed, out, seen, todo, deleted, lost
vars == <<pc, dest, part, holds, failed, out, seen, todo, deleted, lost>>

None == "none"
Procs == Getters \cup Clearers
Names == Getters
\* staged_path: mkstemp gives every writer its own sibling name.
Name(p) == p
NoPart == [c |-> "none", own |-> None]
Absent == [c |-> "absent", gen |-> None]
Fails(p, stage) == p \in FailProcs /\ FailStage = stage

\* clear_cache: what "*.csv" and "annotations_*.zip" match.
Matches == IF dest.c = "absent" THEN {} ELSE {"dest"}

Init ==
  /\ pc = [p \in Procs |-> IF p \in Getters THEN "load" ELSE "glob"]
  /\ dest = [c |-> InitDest, gen |-> None]
  /\ part = [n \in Names |-> NoPart]
  /\ holds = [p \in Getters |-> FALSE]
  /\ failed = [p \in Getters |-> FALSE]
  /\ out = [p \in Getters |-> "none"]
  /\ seen = [p \in Getters |-> "none"]
  /\ todo = [c \in Clearers |-> {}]
  /\ deleted = [c \in Clearers |-> 0]
  /\ lost = FALSE

GLoad(p) ==
  /\ pc[p] = "load"
  /\ seen' = [seen EXCEPT ![p] = dest.c]
  /\ IF dest.c = "good"
       THEN /\ out' = [out EXCEPT ![p] = "hit"]
            /\ pc' = [pc EXCEPT ![p] = "done"]
       ELSE /\ pc' = [pc EXCEPT ![p] = "fetch"]
            /\ UNCHANGED out
  /\ UNCHANGED <<dest, part, holds, failed, todo, deleted, lost>>

GFetch(p) ==
  /\ pc[p] = "fetch"
  /\ IF Fails(p, "fetch")
       THEN /\ out' = [out EXCEPT ![p] = "fetcherr"]
            /\ pc' = [pc EXCEPT ![p] = "done"]
       ELSE /\ pc' = [pc EXCEPT ![p] = "mk"]
            /\ UNCHANGED out
  /\ UNCHANGED <<dest, part, holds, failed, seen, todo, deleted, lost>>

GMk(p) ==
  /\ pc[p] = "mk"
  /\ IF Fails(p, "create")
       THEN /\ out' = [out EXCEPT ![p] = "warned"]
            /\ pc' = [pc EXCEPT ![p] = "done"]
            /\ UNCHANGED <<part, holds>>
       ELSE /\ part' = [part EXCEPT ![Name(p)] = [c |-> "empty", own |-> p]]
            /\ holds' = [holds EXCEPT ![p] = TRUE]
            /\ pc' = [pc EXCEPT ![p] = "open"]
            /\ UNCHANGED out
  /\ UNCHANGED <<dest, failed, seen, todo, deleted, lost>>

\* ZipFile(name, "w") truncates the file, or creates it if the name is gone.
GOpen(p) ==
  /\ pc[p] = "open"
  /\ part' = [part EXCEPT ![Name(p)] = [c |-> "partial", own |-> p]]
  /\ pc' = [pc EXCEPT ![p] = "write"]
  /\ UNCHANGED <<dest, holds, failed, out, seen, todo, deleted, lost>>

\* Writes go through the open descriptor: they cannot recreate a name that
\* was unlinked meanwhile.
GWrite(p) ==
  /\ pc[p] = "write"
  /\ IF Fails(p, "write")
       THEN /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ pc' = [pc EXCEPT ![p] = "cleanup"]
            /\ UNCHANGED part
       ELSE /\ part' = IF part[Name(p)].c = "none" THEN part
                       ELSE [part EXCEPT ![Name(p)] = [c |-> "good", own |-> p]]
            /\ pc' = [pc EXCEPT ![p] = "replace"]
            /\ UNCHANGED failed
  /\ UNCHANGED <<dest, holds, out, seen, todo, deleted, lost>>

GReplace(p) ==
  /\ pc[p] = "replace"
  /\ pc' = [pc EXCEPT ![p] = "cleanup"]
  /\ IF Fails(p, "replace")
       THEN /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ UNCHANGED <<dest, part, holds, lost>>
       ELSE IF part[Name(p)].c = "none"
         THEN \* FileNotFoundError from os.replace: the sibling vanished.
              /\ failed' = [failed EXCEPT ![p] = TRUE]
              /\ lost' = TRUE
              /\ UNCHANGED <<dest, part, holds>>
         ELSE /\ dest' = [c |-> part[Name(p)].c, gen |-> p]
              /\ part' = [part EXCEPT ![Name(p)] = NoPart]
              /\ holds' = [holds EXCEPT ![p] = FALSE]
              /\ UNCHANGED <<failed, lost>>
  /\ UNCHANGED <<out, seen, todo, deleted>>

GCleanup(p) ==
  /\ pc[p] = "cleanup"
  /\ part' = [part EXCEPT ![Name(p)] = NoPart]
  /\ holds' = [holds EXCEPT ![p] = FALSE]
  /\ out' = [out EXCEPT ![p] = IF failed[p] THEN "warned" ELSE "saved"]
  /\ failed' = [failed EXCEPT ![p] = FALSE]
  /\ pc' = [pc EXCEPT ![p] = "done"]
  /\ UNCHANGED <<dest, seen, todo, deleted, lost>>

CGlob(c) ==
  /\ pc[c] = "glob"
  /\ todo' = [todo EXCEPT ![c] = Matches]
  /\ pc' = [pc EXCEPT ![c] = IF Matches = {} THEN "done" ELSE "unlink"]
  /\ UNCHANGED <<dest, part, holds, failed, out, seen, deleted, lost>>

\* A name that is gone raises FileNotFoundError, which clear_cache logs.
CUnlink(c) ==
  /\ pc[c] = "unlink"
  /\ \E f \in todo[c] :
       /\ todo' = [todo EXCEPT ![c] = @ \ {f}]
       /\ pc' = [pc EXCEPT ![c] = IF todo[c] = {f} THEN "done" ELSE "unlink"]
       /\ IF f = "dest"
            THEN /\ dest' = Absent
                 /\ deleted' = [deleted EXCEPT ![c] =
                                  IF dest.c = "absent" THEN @ ELSE @ + 1]
                 /\ UNCHANGED part
            ELSE /\ part' = [part EXCEPT ![f] = NoPart]
                 /\ deleted' = [deleted EXCEPT ![c] =
                                  IF part[f].c = "none" THEN @ ELSE @ + 1]
                 /\ UNCHANGED dest
  /\ UNCHANGED <<holds, failed, out, seen, lost>>

GStep(p) == \/ GLoad(p) \/ GFetch(p) \/ GMk(p) \/ GOpen(p) \/ GWrite(p)
            \/ GReplace(p) \/ GCleanup(p)
CStep(c) == CGlob(c) \/ CUnlink(c)

AllDone == \A p \in Procs : pc[p] = "done"

\* Stutter at the end so that termination is not reported as a deadlock.
Terminated == AllDone /\ UNCHANGED vars

Next == \/ \E p \in Getters : GStep(p)
        \/ \E c \in Clearers : CStep(c)
        \/ Terminated

\* Weak fairness on every process step. No strong fairness is assumed.
Spec == /\ Init /\ [][Next]_vars
        /\ \A p \in Getters : WF_vars(GStep(p))
        /\ \A c \in Clearers : WF_vars(CStep(c))

--------------------------------------------------------------------------
TypeOK ==
  /\ pc \in [Procs -> {"load", "fetch", "mk", "open", "write", "replace",
                       "cleanup", "glob", "unlink", "done"}]
  /\ dest \in [c : {"absent", "good", "corrupt", "partial"}, gen : Getters \cup {None}]
  /\ part \in [Names -> [c : {"none", "empty", "partial", "good"},
                         own : Getters \cup {None}]]
  /\ holds \in [Getters -> BOOLEAN]
  /\ failed \in [Getters -> BOOLEAN]
  /\ out \in [Getters -> {"none", "hit", "saved", "warned", "fetcherr"}]
  /\ seen \in [Getters -> {"none", "absent", "good", "corrupt", "partial"}]
  /\ todo \in [Clearers -> SUBSET ({"dest"} \cup Names)]
  \* One entry matches the globs, so one clear deletes at most one file.
  /\ deleted \in [Clearers -> 0..1]
  /\ lost \in BOOLEAN

\* Every sibling file belongs to exactly one live staged_path frame, and
\* every such frame still has its file. clear_cache never takes one.
PartOwnership ==
  /\ \A n \in Names : part[n].c # "none" =>
        /\ part[n].own \in Getters
        /\ Name(part[n].own) = n
        /\ holds[part[n].own]
  /\ \A p \in Getters : holds[p] =>
        /\ pc[p] \in {"open", "write", "replace", "cleanup"}
        /\ part[Name(p)].c # "none"
        /\ part[Name(p)].own = p

NoLeak == AllDone => \A n \in Names : part[n] = NoPart

\* The entry never holds a half-written archive, and no reader opens one.
DestNeverPartial == dest.c # "partial"
ReaderSeesWhole == \A p \in Getters : seen[p] # "partial"

\* A save is dropped only when that save's own write raised OSError.
SaveLostOnlyOnFault ==
  /\ ~lost
  /\ \A p \in Getters : out[p] = "warned" =>
        p \in FailProcs /\ FailStage \in {"create", "write", "replace"}

\* At quiescence the entry is absent or complete, unless a corrupt entry
\* predates every process.
QuiescentClean ==
  AllDone => (dest.c \in {"absent", "good"} \/ dest.gen = None)

Termination == <>AllDone
====
