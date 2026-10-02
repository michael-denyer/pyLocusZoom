---------------------------- MODULE RecombPublish ----------------------------
(***************************************************************************)
(* Concurrent publication of the recombination map cache.                   *)
(*                                                                         *)
(* Source: src/pylocuszoom/recombination.py                                *)
(*   _has_complete_maps                    l.147-152                       *)
(*   _holds_only_maps                      l.155-161                       *)
(*   _publish_map_generation               l.164-206                       *)
(*   _stage_archive                        l.209-257                       *)
(*   download_recombination_maps           l.260-301                       *)
(*   download_canine_recombination_maps    l.304-341                       *)
(*   load_recombination_map                l.344-413                       *)
(*   get_recombination_rate_for_region     l.416-512                       *)
(*   ensure_recomb_maps                    l.515-550                       *)
(* and src/pylocuszoom/_http.py download_file l.152-192, staged_path       *)
(* l.71-87, whose .part file lives inside the writer's private temporary   *)
(* directory and is covered by that directory's cleanup.                   *)
(*                                                                         *)
(* A process is a separate OS process or notebook kernel sharing one cache *)
(* directory. There is no mutex and no condition variable, so there is no  *)
(* lock variable, no wait set and no lost wakeup to model; a process is    *)
(* never blocked, and TLC's deadlock check together with Termination       *)
(* covers "every call returns". The shared state is the filesystem. Each   *)
(* filesystem call is one atomic step under POSIX rename semantics, and a  *)
(* check (exists, is_dir, is_symlink, glob, iterdir) and the call it       *)
(* guards are separate steps, as they are in the code.                     *)
(*                                                                         *)
(* A writer runs one of four entry points, chosen by mode:                 *)
(*   "ensure"        ensure_recomb_maps(), or                              *)
(*                   download_canine_recombination_maps()                  *)
(*   "force"         download_canine_recombination_maps(force=True)        *)
(*   "checked"       download_canine_recombination_maps(output_dir=p)      *)
(*   "forcechecked"  the same with force=True                              *)
(* The two checked modes run _holds_only_maps first. A writer in Readers   *)
(* is get_recombination_rate_for_region: it runs ensure_recomb_maps and    *)
(* then load_recombination_map's exists() and open.                        *)
(*                                                                         *)
(* Each writer owns one tempfile.TemporaryDirectory (tmp) holding a        *)
(* staging directory (stgAt, stg). The staging directory is renamed into   *)
(* place or its files are moved out one by one; the temporary directory is *)
(* removed on every exit path.                                             *)
(*                                                                         *)
(* Failure injection: a writer in FailAt has its download or extraction    *)
(* raise; a writer in ShortAt extracts an incomplete set. Sticky makes     *)
(* unlinking the legacy symlink fail (a sticky cache parent whose symlink  *)
(* another user owns). ExtSteps bounds an external actor that deletes map  *)
(* files and removes the emptied directory (a user clearing the cache).    *)
(***************************************************************************)
EXTENDS Naturals, FiniteSets

CONSTANTS
    Writers,    \* processes
    Readers,    \* the writers that go on to load a map
    Names,      \* the map set, source.filenames; naturals, so sorted() is Min
    Stray,      \* a chr*_recomb.tsv outside the map set
    Other,      \* an entry the chr*_recomb.tsv glob does not match
    InitKinds,  \* what the cache path may be at the start
    Modes,      \* entry points a writer outside Readers may run
    FailAt,     \* writers whose download or extraction raises
    ShortAt,    \* writers whose archive holds an incomplete map set
    Sticky,     \* TRUE: unlinking the legacy symlink is denied
    ExtSteps    \* steps an external actor may take

ASSUME /\ Readers \subseteq Writers /\ FailAt \subseteq Writers
       /\ ShortAt \subseteq Writers
       /\ Names # {} /\ Names \subseteq Nat /\ IsFiniteSet(Names)
       /\ Sticky \in BOOLEAN /\ ExtSteps \in Nat

Globbed == Names \cup {Stray}
Entries == Globbed \cup {Other}
Kinds == {"absent", "file", "dir", "symlink", "dangling"}
AllModes == {"ensure", "force", "checked", "forcechecked"}
Min(S) == CHOOSE x \in S : \A y \in S : x <= y

VARIABLES
    out,       \* the cache path: [kind, files]; a symlink's files are its target's
    mode,      \* entry point each writer runs
    pc,        \* program counter
    outcome,   \* what the writer returns once its temporary directory is gone
    tmp,       \* the writer's TemporaryDirectory: "none", "live" or "cleaned"
    stgAt,     \* its staging directory: "none", "tmp" (inside tmp) or "out"
    stg,       \* map files still inside its staging directory
    todo,      \* names the fallback loop has still to os.replace
    seen,      \* entries a listing returned (iterdir at l.160, glob at l.203)
    rname,     \* the map a reader loads
    swapped,   \* history: a writer has unlinked the legacy symlink
    foreign0,  \* history: the path started as a file or held a non-map entry
    ext        \* steps the external actor has taken

vars == <<out, mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname, swapped,
          foreign0, ext>>

Force(w) == mode[w] \in {"force", "forcechecked"}
Checked(w) == mode[w] \in {"checked", "forcechecked"}

\* What a path lookup that follows symlinks sees.
IsDir == out.kind \in {"dir", "symlink"}
Exists == out.kind \in {"file", "dir", "symlink"}
Visible == IF IsDir THEN out.files ELSE {}
\* _has_complete_maps (l.151-152): the globbed names are exactly the map set.
Complete == (Visible \cap Globbed) = Names

Init ==
    /\ out \in [kind : InitKinds \cap {"absent", "file", "dangling"}, files : {{}}]
               \cup [kind : InitKinds \cap {"dir", "symlink"}, files : SUBSET Entries]
    /\ mode \in [Writers -> Modes \cup {"ensure"}]
    /\ \A w \in Readers : mode[w] = "ensure"
    /\ \A w \in Writers \ Readers : mode[w] \in Modes
    /\ pc = [w \in Writers |-> "start"]
    /\ outcome = [w \in Writers |-> "none"]
    /\ tmp = [w \in Writers |-> "none"]
    /\ stgAt = [w \in Writers |-> "none"]
    /\ stg = [w \in Writers |-> {}]
    /\ todo = [w \in Writers |-> {}]
    /\ seen = [w \in Writers |-> {}]
    /\ rname \in [Writers -> Names]
    /\ \A w \in Writers \ Readers : rname[w] = Min(Names)
    /\ swapped = FALSE
    /\ foreign0 = (out.kind = "file" \/ out.files \cap {Stray, Other} # {})
    /\ ext = 0

Goto(w, label) == pc' = [pc EXCEPT ![w] = label]
\* Leave the with block at l.285: the exit runs before the caller sees o.
Leave(w, o) == /\ Goto(w, "cleanup")
               /\ outcome' = [outcome EXCEPT ![w] = o]

AfterMiss(w) == IF Checked(w) THEN "hexists" ELSE "mktmp"
AfterHit(w) == IF w \in Readers THEN "rexists" ELSE "hit"

\* l.333 `not force and ...` / l.546.
Start(w) ==
    /\ pc[w] = "start"
    /\ Goto(w, IF Force(w) THEN AfterMiss(w) ELSE "cexists")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* _has_complete_maps l.149: path.exists().
CExists(w) ==
    /\ pc[w] = "cexists"
    /\ Goto(w, IF Exists THEN "cglob" ELSE AfterMiss(w))
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* _has_complete_maps l.151-152: the glob and the set comparison.
CGlob(w) ==
    /\ pc[w] = "cglob"
    /\ Goto(w, IF Complete THEN AfterHit(w) ELSE AfterMiss(w))
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* _holds_only_maps l.157: path.exists().
HExists(w) ==
    /\ pc[w] = "hexists"
    /\ Goto(w, IF Exists THEN "hisdir" ELSE "mktmp")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* _holds_only_maps l.159: path.is_dir(); False raises ValidationError (l.336).
HIsDir(w) ==
    /\ pc[w] = "hisdir"
    /\ Goto(w, IF IsDir THEN "hlist" ELSE "invalid")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* _holds_only_maps l.160: path.iterdir(). On a path that is no longer a
\* directory it raises FileNotFoundError or NotADirectoryError, which nothing
\* between here and the caller catches.
HList(w) ==
    /\ pc[w] = "hlist"
    /\ IF IsDir
         THEN /\ Goto(w, "hstat")
              /\ seen' = [seen EXCEPT ![w] = out.files]
         ELSE /\ Goto(w, "crashed")
              /\ UNCHANGED seen
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, rname, swapped,
                   foreign0, ext>>

\* _holds_only_maps l.160: entry.is_file() and the name test for the listed
\* entries. An entry that has vanished since the listing is not a file.
HStat(w) ==
    /\ pc[w] = "hstat"
    /\ Goto(w, IF seen[w] \subseteq Names /\ seen[w] \subseteq Visible
                 THEN "mktmp" ELSE "invalid")
    /\ seen' = [seen EXCEPT ![w] = {}]
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, rname, swapped,
                   foreign0, ext>>

\* download_recombination_maps l.284-285: TemporaryDirectory(dir=parent).
MkTmp(w) ==
    /\ pc[w] = "mktmp"
    /\ tmp' = [tmp EXCEPT ![w] = "live"]
    /\ Goto(w, "download")
    /\ UNCHANGED <<out, mode, outcome, stgAt, stg, todo, seen, rname, swapped,
                   foreign0, ext>>

\* l.288 download_file: raises DataDownloadError for a writer in FailAt,
\* which may instead fail later, inside _stage_archive.
Download(w) ==
    /\ pc[w] = "download"
    /\ \/ /\ Goto(w, "stage")
          /\ UNCHANGED outcome
       \/ /\ w \in FailAt
          /\ Leave(w, "dlfailed")
    /\ UNCHANGED <<out, mode, tmp, stgAt, stg, todo, seen, rname, swapped,
                   foreign0, ext>>

\* l.291-293: staging.mkdir() and _stage_archive. A FailAt writer raises
\* part-way, leaving any subset of the maps in its staging directory.
Stage(w) ==
    /\ pc[w] = "stage"
    /\ stgAt' = [stgAt EXCEPT ![w] = "tmp"]
    /\ IF w \in FailAt
         THEN /\ \E s \in SUBSET Names : stg' = [stg EXCEPT ![w] = s]
              /\ Leave(w, "dlfailed")
         ELSE /\ IF w \in ShortAt
                   THEN \E s \in (SUBSET Names) \ {Names} :
                            stg' = [stg EXCEPT ![w] = s]
                   ELSE stg' = [stg EXCEPT ![w] = Names]
              /\ Goto(w, "pubcheck")
              /\ UNCHANGED outcome
    /\ UNCHANGED <<out, mode, tmp, todo, seen, rname, swapped, foreign0, ext>>

\* _publish_map_generation l.180: _has_complete_maps(staging_dir). The
\* staging directory is private, so its exists() and glob are one step.
PubCheck(w) ==
    /\ pc[w] = "pubcheck"
    /\ IF stg[w] = Names
         THEN Goto(w, "symcheck") /\ UNCHANGED outcome
         ELSE Leave(w, "dlfailed")
    /\ UNCHANGED <<out, mode, tmp, stgAt, stg, todo, seen, rname, swapped,
                   foreign0, ext>>

\* l.187: if output_path.is_symlink():
SymCheck(w) ==
    /\ pc[w] = "symcheck"
    /\ Goto(w, IF out.kind \in {"symlink", "dangling"} THEN "unlink" ELSE "rename")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* l.190 output_path.unlink(): FileNotFoundError once the link is gone,
\* EPERM or EISDIR on a directory, EPERM on a link Sticky protects.
Unlink(w) ==
    /\ pc[w] = "unlink"
    /\ IF out.kind = "file" \/ (out.kind \in {"symlink", "dangling"} /\ ~Sticky)
         THEN /\ out' = [kind |-> "absent", files |-> {}]
              /\ swapped' = (swapped \/ out.kind # "file")
              /\ Goto(w, "rename")
         ELSE /\ Goto(w, "recheck")
              /\ UNCHANGED <<out, swapped>>
    /\ UNCHANGED <<mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   foreign0, ext>>

\* l.193-194: except OSError: if output_path.is_symlink(): raise
Recheck(w) ==
    /\ pc[w] = "recheck"
    /\ IF out.kind \in {"symlink", "dangling"}
         THEN Leave(w, "failed")
         ELSE Goto(w, "rename") /\ UNCHANGED outcome
    /\ UNCHANGED <<out, mode, tmp, stgAt, stg, todo, seen, rname, swapped,
                   foreign0, ext>>

\* l.196 os.rename(staging_dir, output_path): succeeds onto nothing or an
\* empty directory, and the staging directory becomes the cache. It fails
\* onto a non-empty directory, a file or a symlink.
Rename(w) ==
    /\ pc[w] = "rename"
    /\ IF out.kind = "absent" \/ (out.kind = "dir" /\ out.files = {})
         THEN /\ out' = [kind |-> "dir", files |-> stg[w]]
              /\ stg' = [stg EXCEPT ![w] = {}]
              /\ stgAt' = [stgAt EXCEPT ![w] = "out"]
              /\ Leave(w, "done")
         ELSE /\ Goto(w, "isdir")
              /\ UNCHANGED <<out, stg, stgAt, outcome>>
    /\ UNCHANGED <<mode, tmp, todo, seen, rname, swapped, foreign0, ext>>

\* l.199-200: except OSError: if not output_path.is_dir(): raise
IsDirCheck(w) ==
    /\ pc[w] = "isdir"
    /\ IF IsDir
         THEN /\ Goto(w, "replace")
              /\ todo' = [todo EXCEPT ![w] = Names]
              /\ UNCHANGED outcome
         ELSE /\ Leave(w, "failed")
              /\ UNCHANGED todo
    /\ UNCHANGED <<out, mode, tmp, stgAt, stg, seen, rname, swapped,
                   foreign0, ext>>

\* l.201-202: os.replace(staging_dir / name, output_path / name) in sorted
\* order, one file per step. It raises when the source is missing or the
\* target directory is gone.
Replace(w) ==
    /\ pc[w] = "replace"
    /\ IF todo[w] = {}
         THEN /\ Goto(w, "sglob")
              /\ UNCHANGED <<out, stg, todo, outcome>>
         ELSE LET n == Min(todo[w]) IN
              IF IsDir /\ n \in stg[w]
                THEN /\ out' = [out EXCEPT !.files = @ \cup {n}]
                     /\ stg' = [stg EXCEPT ![w] = @ \ {n}]
                     /\ todo' = [todo EXCEPT ![w] = @ \ {n}]
                     /\ UNCHANGED <<pc, outcome>>
                ELSE /\ Leave(w, "failed")
                     /\ UNCHANGED <<out, stg, todo>>
    /\ UNCHANGED <<mode, tmp, stgAt, seen, rname, swapped, foreign0, ext>>

\* l.203-204: output_path.glob("chr*_recomb.tsv") filtered to non-map names.
\* A glob of a missing directory yields nothing.
StrayGlob(w) ==
    /\ pc[w] = "sglob"
    /\ IF Stray \in Visible
         THEN /\ seen' = [seen EXCEPT ![w] = {Stray}]
              /\ Goto(w, "sunlink")
              /\ UNCHANGED outcome
         ELSE /\ Leave(w, "done")
              /\ UNCHANGED seen
    /\ UNCHANGED <<out, mode, tmp, stgAt, stg, todo, rname, swapped,
                   foreign0, ext>>

\* l.205: stray.unlink(missing_ok=True).
StrayUnlink(w) ==
    /\ pc[w] = "sunlink"
    /\ out' = IF IsDir THEN [out EXCEPT !.files = @ \ seen[w]] ELSE out
    /\ seen' = [seen EXCEPT ![w] = {}]
    /\ Leave(w, "done")
    /\ UNCHANGED <<mode, tmp, stgAt, stg, todo, rname, swapped, foreign0, ext>>

\* TemporaryDirectory.__exit__ at l.285, on return and on every exception:
\* removes the archive and whatever is left of the staging directory. After
\* a successful rename the staging directory is no longer inside it.
Cleanup(w) ==
    /\ pc[w] = "cleanup"
    /\ tmp' = [tmp EXCEPT ![w] = "cleaned"]
    /\ stg' = [stg EXCEPT ![w] = {}]
    /\ Goto(w, IF outcome[w] = "done" THEN AfterHit(w) ELSE outcome[w])
    /\ UNCHANGED <<out, mode, outcome, stgAt, todo, seen, rname, swapped,
                   foreign0, ext>>

\* load_recombination_map l.370: map_file.exists().
ReadExists(w) ==
    /\ pc[w] = "rexists"
    /\ Goto(w, IF rname[w] \in Visible THEN "ropen" ELSE "notfound")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

\* load_recombination_map l.383: pd.read_csv opens the file.
ReadOpen(w) ==
    /\ pc[w] = "ropen"
    /\ Goto(w, IF rname[w] \in Visible THEN "read" ELSE "unreadable")
    /\ UNCHANGED <<out, mode, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0, ext>>

Proc(w) ==
    \/ Start(w) \/ CExists(w) \/ CGlob(w)
    \/ HExists(w) \/ HIsDir(w) \/ HList(w) \/ HStat(w)
    \/ MkTmp(w) \/ Download(w) \/ Stage(w) \/ PubCheck(w)
    \/ SymCheck(w) \/ Unlink(w) \/ Recheck(w) \/ Rename(w) \/ IsDirCheck(w)
    \/ Replace(w) \/ StrayGlob(w) \/ StrayUnlink(w) \/ Cleanup(w)
    \/ ReadExists(w) \/ ReadOpen(w)

\* An external actor clearing the cache: it deletes one entry, or removes
\* the directory once it is empty. It has no fairness requirement.
External ==
    /\ ext < ExtSteps
    /\ out.kind = "dir"
    /\ \/ \E n \in out.files : out' = [out EXCEPT !.files = @ \ {n}]
       \/ out.files = {} /\ out' = [kind |-> "absent", files |-> {}]
    /\ ext' = ext + 1
    /\ UNCHANGED <<mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname,
                   swapped, foreign0>>

TerminalPcs == {"hit", "done", "failed", "dlfailed", "invalid", "crashed",
                "read", "notfound", "unreadable"}
AllTerminal == \A w \in Writers : pc[w] \in TerminalPcs

\* Stutter at the end so that termination is not reported as a deadlock.
Terminated == AllTerminal /\ UNCHANGED vars

Next == (\E w \in Writers : Proc(w)) \/ External \/ Terminated

\* Weak fairness on every process step. No step waits on another process,
\* so no strong fairness assumption is needed.
Spec == Init /\ [][Next]_vars /\ \A w \in Writers : WF_vars(Proc(w))

-----------------------------------------------------------------------------
TypeOK ==
    /\ out \in [kind : Kinds, files : SUBSET Entries]
    /\ out.kind \in {"absent", "file", "dangling"} => out.files = {}
    /\ mode \in [Writers -> AllModes]
    /\ pc \in [Writers -> TerminalPcs \cup
                 {"start", "cexists", "cglob", "hexists", "hisdir", "hlist",
                  "hstat", "mktmp", "download", "stage", "pubcheck",
                  "symcheck", "unlink", "recheck", "rename", "isdir",
                  "replace", "sglob", "sunlink", "cleanup", "rexists",
                  "ropen"}]
    /\ outcome \in [Writers -> {"none", "done", "failed", "dlfailed"}]
    /\ tmp \in [Writers -> {"none", "live", "cleaned"}]
    /\ stgAt \in [Writers -> {"none", "tmp", "out"}]
    /\ stg \in [Writers -> SUBSET Names]
    /\ todo \in [Writers -> SUBSET Names]
    /\ seen \in [Writers -> SUBSET Entries]
    /\ rname \in [Writers -> Names]
    /\ swapped \in BOOLEAN /\ foreign0 \in BOOLEAN
    /\ ext \in 0..ExtSteps

\* The call returned a usable cache: a hit, a publish, or any later reader
\* step.
Succeeded(w) == pc[w] \in {"hit", "done", "rexists", "ropen", "read",
                           "notfound", "unreadable"}
Failing(w) == pc[w] = "failed" \/ (pc[w] = "cleanup" /\ outcome[w] = "failed")

\* Resource ownership. Staged files exist only inside a live temporary
\* directory; a staging directory that became the cache is never used as a
\* source again; the fallback loop only takes files that are still staged.
StagingOwned ==
    \A w \in Writers :
        /\ stg[w] # {} => tmp[w] = "live" /\ stgAt[w] = "tmp"
        /\ stgAt[w] # "none" => tmp[w] # "none"
        /\ pc[w] \in {"isdir", "replace", "sglob", "sunlink"} => stgAt[w] = "tmp"
        /\ pc[w] = "replace" => todo[w] \subseteq stg[w]
        /\ stgAt[w] = "out" => Succeeded(w) \/ pc[w] = "cleanup"

\* No writer returns, by any path, with its temporary directory still in
\* the cache parent.
NoLeak ==
    \A w \in Writers :
        pc[w] \in TerminalPcs \cup {"rexists", "ropen"} => tmp[w] # "live"

\* At most one staging directory ever becomes the cache directory: no
\* writer's rename discards a directory another writer published.
SinglePublisher == Cardinality({w \in Writers : stgAt[w] = "out"}) <= 1

\* "a writer that loses the race to another still succeeds": publication
\* raises only where no map set can be published, a regular file at the
\* cache path or a legacy symlink the writer may not remove.
WriterSucceeds ==
    \A w \in Writers :
        Failing(w) => \/ out.kind = "file"
                      \/ Sticky /\ out.kind \in {"symlink", "dangling"}

\* A download raises only where a failure was injected.
FailsOnlyWhenInjected ==
    \A w \in Writers :
        pc[w] = "dlfailed" \/ outcome[w] = "dlfailed" => w \in FailAt \cup ShortAt

\* Once every call has returned and one of them succeeded, the cache is a
\* complete map set that the next ensure_recomb_maps accepts as a hit.
Converges == AllTerminal /\ (\E w \in Writers : Succeeded(w)) => Complete

\* A writer never writes through a legacy symlink into a directory the
\* package does not own ("its target is not ours").
NeverWritesThroughSymlink ==
    \A w \in Writers :
        pc[w] \in {"replace", "sglob", "sunlink"} => out.kind # "symlink"

\* A checked writer never runs the fallback over a directory holding
\* anything besides the map set.
CheckedNeverClobbers ==
    \A w \in Writers :
        Checked(w) /\ pc[w] \in {"replace", "sglob", "sunlink"}
            => Visible \cap {Stray, Other} = {}

\* Strict claims. They hold when the cache path does not start as a legacy
\* symlink.

\* "a reader finds the old file or the new one, never a gap"
NoGap == \A w \in Writers : pc[w] \notin {"notfound", "unreadable"}

\* From the first successful return on, the cache stays complete.
StaysComplete == (\E w \in Writers : Succeeded(w)) => Complete

\* Only the documented exceptions leave download_canine_recombination_maps.
NoUndocumentedError == \A w \in Writers : pc[w] # "crashed"

\* ValidationError is raised only for a path that held something else.
ValidationJustified == \A w \in Writers : pc[w] = "invalid" => foreign0

\* The same four claims with the documented exception: each can break only
\* after a writer has unlinked a legacy symlink, and the cache is incomplete
\* only in the window between that unlink and the rename that follows it.
NoGapUnlessSwapped ==
    \A w \in Writers : pc[w] \in {"notfound", "unreadable"} => swapped
StaysCompleteUnlessSwapping ==
    (\E w \in Writers : Succeeded(w))
        => Complete \/ (swapped /\ out.kind = "absent")
NoUndocumentedErrorUnlessSwapped ==
    \A w \in Writers : pc[w] = "crashed" => swapped
ValidationJustifiedUnlessSwapped ==
    \A w \in Writers : pc[w] = "invalid" => foreign0 \/ swapped

\* Every call returns: each writer reaches a terminal pc and each reader
\* finishes its load.
Termination == <>AllTerminal
=============================================================================
