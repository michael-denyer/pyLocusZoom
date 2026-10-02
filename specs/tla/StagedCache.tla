---- MODULE StagedCache ----
(* The liftover-chain cache: a download is parsed while it is still a private
   sibling and replaces the cached file only if it parses, for separate OS
   processes sharing one cache directory. StagedGeneCache.tla models the gene
   cache over the same writer.

   Source modelled (src/pylocuszoom/):
     _http.py::staged_path          mkstemp sibling, yield, os.replace, finally unlink
     _http.py::stream_file          (abstracted, see below)
     _http.py::_with_retries        (abstracted, see below)
     _liftover.py::load_chain       and _liftover.py::_parse_chain
     _liftover.py::chain_lifter     resolve the URL and the cache path
     _liftover.py::_cached_chain    @lru_cache(maxsize=4) per process:
                                    load, else download
     _liftover.py::_download_chain  stage, download, parse, publish

   One action is one filesystem call. A check and the call it guards are
   separate actions, so every interleaving between them is explored.

     pc        code                                             action
     load      _cached_chain: memo, load_chain(path)            Load
     mk        _download_chain: mkdir, staged_path mkstemp      Mk
     open      _download_chain: stream_file starts              Open
     stream    _download_chain: stream_file ends or raises      Stream
     validate  _download_chain: _parse_chain(partial)           Validate
     replace   staged_path: os.replace(partial_path, dest)      Replace
     cleanup   staged_path: partial_path.unlink(missing_ok)     Cleanup

   No action unlinks the destination: the only writer of dest is Replace, and
   Validate stands between every download and it.

   load_chain opens the file once and parses from that descriptor, so one
   Load action observes one content. An absent file is FileNotFoundError,
   which _parse_chain turns into ValidationError like any unreadable content.
   Validate parses the sibling the way the destination would be opened
   (gzipped by the destination's suffix); content that does not parse raises,
   staged_path skips os.replace and removes the sibling, and the call raises
   DataDownloadError. A partially written file is modelled as parsing
   "successfully" (a truncated plain chain can parse), so LifterGood and
   DestNeverPartial depend on the staging protocol and not on the parser.

   `memo` is _cached_chain's lru_cache: once a call has returned a lifter in
   a process, later calls return it without touching the file. It is filled
   by a load of the cached file and by a validated download alike; `lifter`
   is the lifter parsed from the sibling and not yet returned.

   Failure injection is by constants. Download n of process p (counted across
   its calls) raises at FailStage when p \in FailProcs and n \in FailDl; it
   yields content that does not parse when p \in CorruptProcs and
   n \in CorruptDl. _with_retries is abstracted: a retry reopens a private
   sibling with "wb", which no other process can observe, so only the final
   outcome of the retry loop (complete or raise) is a step.

   Damage models what the code does not do: with Damage # "none" the
   environment may, once, at any moment, remove the cached file ("absent")
   or replace it with unparseable content ("corrupt"), as a user clearing
   the cache, disk damage or an older version of the package would.

   Ownership: part[n] is the sibling file named n and records which process
   created it; holds[p] says p's live staged_path frame owns a sibling.
   PartOwnership is the partition, NoLeak the terminal condition.

   There is no mutex and no condition variable in this protocol, so the
   model has no lock variable, no wait set and no spurious wakeups, and a
   lost wakeup cannot occur. Every action is always enabled once its pc is
   reached; liveness (Termination) therefore checks that every call ends,
   under weak fairness on each process's steps. Damage is not fair: it may
   never happen. *)
EXTENDS Naturals

CONSTANTS Procs, Calls, InitDest, FailProcs, FailDl, FailStage,
          CorruptProcs, CorruptDl, Damage

ASSUME /\ "none" \notin Procs
       /\ Calls \in 1..2
       /\ InitDest \in {"absent", "good", "corrupt"}
       /\ FailProcs \subseteq Procs /\ CorruptProcs \subseteq Procs
       /\ FailStage \in {"none", "create", "stream", "replace"}
       /\ Damage \in {"none", "absent", "corrupt"}

VARIABLES pc, failed, dl, calls, dest, part, holds, lifter, memo, res, blame,
          clobber, lost, damaged
vars == <<pc, failed, dl, calls, dest, part, holds, lifter, memo, res, blame,
          clobber, lost, damaged>>

None == "none"
\* A call downloads at most once.
MaxDl == Calls
Names == Procs
\* staged_path: mkstemp gives every writer its own sibling name.
Name(p) == p
NoPart == [c |-> "none", own |-> None]

Fails(p, n, stage) == p \in FailProcs /\ n \in FailDl /\ FailStage = stage
Corrupts(p, n) == p \in CorruptProcs /\ n \in CorruptDl

Init ==
  /\ pc = [p \in Procs |-> "idle"]
  /\ failed = [p \in Procs |-> FALSE]
  /\ dl = [p \in Procs |-> 0]
  /\ calls = [p \in Procs |-> 0]
  /\ dest = [c |-> InitDest, gen |-> None]
  /\ part = [n \in Names |-> NoPart]
  /\ holds = [p \in Procs |-> FALSE]
  /\ lifter = [p \in Procs |-> "none"]
  /\ memo = [p \in Procs |-> "none"]
  /\ res = [p \in Procs |-> "none"]
  /\ blame = [p \in Procs |-> FALSE]
  /\ clobber = FALSE
  /\ lost = FALSE
  /\ damaged = FALSE

\* The call returns or raises: back to the caller.
Finish(p, r) ==
  /\ res' = [res EXCEPT ![p] = r]
  /\ calls' = [calls EXCEPT ![p] = @ + 1]
  /\ pc' = [pc EXCEPT ![p] = "idle"]

Start(p) ==
  /\ pc[p] = "idle" /\ calls[p] < Calls
  /\ pc' = [pc EXCEPT ![p] = "load"]
  /\ res' = [res EXCEPT ![p] = "none"]
  /\ blame' = [blame EXCEPT ![p] = FALSE]
  /\ lifter' = [lifter EXCEPT ![p] = "none"]
  /\ UNCHANGED <<failed, dl, calls, dest, part, holds, memo, clobber, lost,
                 damaged>>

\* _cached_chain returns the memoised lifter; otherwise load_chain parses the cached file,
\* and anything that does not parse sends the call to the download.
Load(p) ==
  /\ pc[p] = "load"
  /\ IF memo[p] # "none"
       THEN Finish(p, "ok") /\ UNCHANGED memo
       ELSE IF dest.c \in {"good", "partial"}
         THEN /\ memo' = [memo EXCEPT ![p] = dest.c]
              /\ Finish(p, "ok")
         ELSE /\ pc' = [pc EXCEPT ![p] = "mk"]
              /\ UNCHANGED <<memo, res, calls>>
  /\ UNCHANGED <<failed, dl, dest, part, holds, lifter, blame, clobber, lost,
                 damaged>>

Mk(p) ==
  /\ pc[p] = "mk"
  /\ dl' = [dl EXCEPT ![p] = @ + 1]
  /\ IF Fails(p, dl[p] + 1, "create")
       THEN /\ Finish(p, "dde")
            /\ blame' = [blame EXCEPT ![p] = TRUE]
            /\ UNCHANGED <<part, holds>>
       ELSE /\ part' = [part EXCEPT ![Name(p)] = [c |-> "empty", own |-> p]]
            /\ holds' = [holds EXCEPT ![p] = TRUE]
            /\ pc' = [pc EXCEPT ![p] = "open"]
            /\ UNCHANGED <<res, calls, blame>>
  /\ UNCHANGED <<failed, dest, lifter, memo, clobber, lost, damaged>>

Open(p) ==
  /\ pc[p] = "open"
  /\ part' = [part EXCEPT ![Name(p)] = [c |-> "partial", own |-> p]]
  /\ pc' = [pc EXCEPT ![p] = "stream"]
  /\ UNCHANGED <<failed, dl, calls, dest, holds, lifter, memo, res, blame,
                 clobber, lost, damaged>>

Stream(p) ==
  /\ pc[p] = "stream"
  /\ IF Fails(p, dl[p], "stream")
       THEN /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ blame' = [blame EXCEPT ![p] = TRUE]
            /\ pc' = [pc EXCEPT ![p] = "cleanup"]
            /\ UNCHANGED part
       ELSE /\ part' = IF part[Name(p)].c = "none" THEN part
                       ELSE [part EXCEPT ![Name(p)] =
                              [c |-> IF Corrupts(p, dl[p]) THEN "corrupt" ELSE "good",
                               own |-> p]]
            /\ blame' = [blame EXCEPT ![p] = @ \/ Corrupts(p, dl[p])]
            /\ pc' = [pc EXCEPT ![p] = "validate"]
            /\ UNCHANGED failed
  /\ UNCHANGED <<dl, calls, dest, holds, lifter, memo, res, clobber, lost,
                 damaged>>

\* _parse_chain(partial) parses the sibling. Content that does not parse, or a sibling that
\* is gone, raises inside the staged_path block: no os.replace, then cleanup.
Validate(p) ==
  /\ pc[p] = "validate"
  /\ IF part[Name(p)].c \in {"good", "partial"}
       THEN /\ lifter' = [lifter EXCEPT ![p] = part[Name(p)].c]
            /\ pc' = [pc EXCEPT ![p] = "replace"]
            /\ UNCHANGED failed
       ELSE /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ pc' = [pc EXCEPT ![p] = "cleanup"]
            /\ UNCHANGED lifter
  /\ UNCHANGED <<dl, calls, dest, part, holds, memo, res, blame, clobber,
                 lost, damaged>>

Replace(p) ==
  /\ pc[p] = "replace"
  /\ pc' = [pc EXCEPT ![p] = "cleanup"]
  /\ IF Fails(p, dl[p], "replace")
       THEN /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ blame' = [blame EXCEPT ![p] = TRUE]
            /\ UNCHANGED <<dest, part, holds, lost, clobber>>
       ELSE IF part[Name(p)].c = "none"
         THEN \* FileNotFoundError from os.replace: the sibling vanished.
              /\ failed' = [failed EXCEPT ![p] = TRUE]
              /\ lost' = TRUE
              /\ UNCHANGED <<dest, part, holds, blame, clobber>>
         ELSE /\ dest' = [c |-> part[Name(p)].c, gen |-> p]
              /\ clobber' = (clobber \/ (dest.c = "good" /\ part[Name(p)].c # "good"))
              /\ part' = [part EXCEPT ![Name(p)] = NoPart]
              /\ holds' = [holds EXCEPT ![p] = FALSE]
              /\ UNCHANGED <<failed, blame, lost>>
  /\ UNCHANGED <<dl, calls, lifter, memo, res, damaged>>

\* The finally clause, then _download_chain returns the lifter or the error
\* propagates.
Cleanup(p) ==
  /\ pc[p] = "cleanup"
  /\ part' = [part EXCEPT ![Name(p)] = NoPart]
  /\ holds' = [holds EXCEPT ![p] = FALSE]
  /\ failed' = [failed EXCEPT ![p] = FALSE]
  /\ IF failed[p]
       THEN Finish(p, "dde") /\ UNCHANGED memo
       ELSE /\ memo' = [memo EXCEPT ![p] = lifter[p]]
            /\ Finish(p, "ok")
  /\ UNCHANGED <<dl, dest, lifter, blame, clobber, lost, damaged>>

Step(p) == \/ Start(p) \/ Load(p) \/ Mk(p) \/ Open(p) \/ Stream(p)
           \/ Validate(p) \/ Replace(p) \/ Cleanup(p)

AllDone == \A p \in Procs : pc[p] = "idle" /\ calls[p] = Calls

\* The environment removes or damages the cached file, at most once.
DamageDest ==
  /\ Damage # "none" /\ ~damaged /\ ~AllDone
  /\ dest' = [c |-> Damage, gen |-> None]
  /\ damaged' = TRUE
  /\ UNCHANGED <<pc, failed, dl, calls, part, holds, lifter, memo, res, blame,
                 clobber, lost>>

\* Stutter at the end so that termination is not reported as a deadlock.
Terminated == AllDone /\ UNCHANGED vars

Next == (\E p \in Procs : Step(p)) \/ DamageDest \/ Terminated

\* Weak fairness on every process step. No strong fairness is assumed.
Spec == Init /\ [][Next]_vars /\ \A p \in Procs : WF_vars(Step(p))

--------------------------------------------------------------------------
TypeOK ==
  /\ pc \in [Procs -> {"idle", "load", "mk", "open", "stream", "validate",
                       "replace", "cleanup"}]
  /\ failed \in [Procs -> BOOLEAN]
  /\ dl \in [Procs -> 0..MaxDl]
  /\ calls \in [Procs -> 0..Calls]
  /\ dest \in [c : {"absent", "good", "corrupt", "partial"}, gen : Procs \cup {None}]
  /\ part \in [Names -> [c : {"none", "empty", "partial", "good", "corrupt"},
                         own : Procs \cup {None}]]
  /\ holds \in [Procs -> BOOLEAN]
  /\ lifter \in [Procs -> {"none", "good", "partial"}]
  /\ memo \in [Procs -> {"none", "good", "partial"}]
  /\ res \in [Procs -> {"none", "ok", "dde"}]
  /\ blame \in [Procs -> BOOLEAN]
  /\ clobber \in BOOLEAN /\ lost \in BOOLEAN /\ damaged \in BOOLEAN

\* Every sibling file belongs to exactly one live staged_path frame, and
\* every such frame still has its file.
PartOwnership ==
  /\ \A n \in Names : part[n].c # "none" =>
        /\ part[n].own \in Procs
        /\ Name(part[n].own) = n
        /\ holds[part[n].own]
  /\ \A p \in Procs : holds[p] =>
        /\ pc[p] \in {"open", "stream", "validate", "replace", "cleanup"}
        /\ part[Name(p)].c # "none"
        /\ part[Name(p)].own = p

\* os.replace always finds the sibling it was given.
PartNeverVanishes == ~lost

NoLeak == AllDone => \A n \in Names : part[n] = NoPart

\* The destination never holds a half-written file.
DestNeverPartial == dest.c # "partial"

\* chain_lifter never hands out a lifter built from anything but a complete,
\* parseable chain.
LifterGood == \A p \in Procs : /\ memo[p] \in {"none", "good"}
                              /\ lifter[p] \in {"none", "good"}

\* No process removes a good chain or replaces it with content that does not
\* parse. (Before the fix an unlink in chain_lifter did the first and every
\* corrupt download did the second.)
NoGoodUnlink == ~clobber

\* A call raises only if its own download failed or was corrupt.
NoSpuriousFailure == \A p \in Procs : res[p] = "dde" => blame[p]

\* A process that already holds a lifter does not raise.
CachedLifterServed == \A p \in Procs : res[p] = "dde" => memo[p] = "none"

\* Whatever a process published parsed when it was published.
PublishedParses == dest.gen # None => dest.c = "good"

\* At quiescence the destination is absent or good, unless the unparseable
\* file came from outside (InitDest or Damage).
QuiescentClean ==
  AllDone => (dest.c \in {"absent", "good"} \/ dest.gen = None)

Termination == <>AllDone
====
