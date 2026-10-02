---- MODULE StagedCache ----
(* The liftover-chain cache: temp-file-then-replace publication and the
   refetch loop around it, for separate OS processes sharing one cache
   directory. StagedGeneCache.tla models the gene cache over the same writer.

   Source modelled (src/pylocuszoom/):
     _http.py:71-87      staged_path: mkstemp sibling, yield, os.replace, finally unlink
     _http.py:152-192    download_file: staged_path around the retried stream
     _http.py:195-214    _stream_to: open(partial, "wb"), write chunks
     _http.py:45-68      _with_retries (abstracted, see below)
     _liftover.py:74-107 load_chain, @lru_cache(maxsize=4) per process and path
     _liftover.py:110-150 chain_lifter: exists, download, load, refetch, load, unlink
     _liftover.py:153-160 _download_chain: OSError becomes DataDownloadError

   One action is one filesystem call. A check and the call it guards are
   separate actions, so every interleaving between them is explored.

     pc        code                                        action
     exists    _liftover.py:134 path.exists()              Exists
     mk        _liftover.py:156, _http.py:79-81 mkstemp    Mk
     open      _http.py:201 open(partial, "wb")            Open
     stream    _http.py:211-213 chunks to end, or raise    Stream
     replace   _http.py:85 os.replace(partial, dest)       Replace
     cleanup   _http.py:87 partial.unlink(missing_ok)      Cleanup
     load      _liftover.py:137 and :142 load_chain        Load
     unlink    _liftover.py:145 path.unlink, :150 raise    Unlink

   `phase` is 1 for the download and load at :135/:137 and 2 for :140/:142.
   load_chain opens the file once and parses from that descriptor, so one
   Load action observes one content. An absent file is FileNotFoundError,
   which :103 turns into ValidationError like any unreadable content. A
   partially written file is modelled as loading "successfully" (a truncated
   plain chain can parse), so LifterGood depends on the staging protocol and
   not on the parser. `memo` is the lru_cache: once a process has a lifter,
   load_chain returns it without touching the file.

   Failure injection is by constants. Download n of process p (counted across
   its calls) raises at FailStage when p \in FailProcs and n \in FailDl; it
   yields content that does not parse when p \in CorruptProcs and
   n \in CorruptDl. _with_retries is abstracted: a retry reopens the private
   sibling with "wb", which no other process can observe, so only the final
   outcome of the retry loop (complete or raise) is a step.

   Ownership: part[n] is the sibling file named n and records which process
   created it; holds[p] says p's live staged_path frame owns a sibling.
   PartOwnership is the partition, NoLeak the terminal condition.

   There is no mutex and no condition variable in this protocol, so the
   model has no lock variable, no wait set and no spurious wakeups, and a
   lost wakeup cannot occur. Every action is always enabled once its pc is
   reached; liveness (Termination) therefore checks that the refetch loop is
   bounded, under weak fairness on each process's steps.

   Witness runs: with Witness # "none" the scheduler follows Sched, a fixed
   interleaving, before running freely. All data choices are constants, so
   the scheduled prefix is the only behaviour and "eventually Goal" passing
   proves the state is reachable. These are the expected failures of
   NoGoodUnlink, NoSpuriousFailure, QuiescentClean and CachedLifterServed
   under a source that returns corrupt content; the matrix checks those four
   claims only where they hold. *)
EXTENDS Naturals, Sequences

CONSTANTS Procs, Calls, InitDest, FailProcs, FailDl, FailStage,
          CorruptProcs, CorruptDl, Witness

ASSUME /\ "none" \notin Procs
       /\ Calls \in 1..2
       /\ InitDest \in {"absent", "good", "corrupt"}
       /\ FailProcs \subseteq Procs /\ CorruptProcs \subseteq Procs
       /\ FailStage \in {"none", "create", "stream", "replace"}

VARIABLES pc, phase, failed, dl, calls, dest, part, holds, memo, res, blame,
          clobber, lost, clock
vars == <<pc, phase, failed, dl, calls, dest, part, holds, memo, res, blame,
          clobber, lost, clock>>

None == "none"
MaxDl == 2 * Calls
Names == Procs
\* _http.py:79-81: mkstemp gives every writer its own sibling name.
Name(p) == p
NoPart == [c |-> "none", own |-> None]
Absent == [c |-> "absent", gen |-> None]

Fails(p, n, stage) == p \in FailProcs /\ n \in FailDl /\ FailStage = stage
Corrupts(p, n) == p \in CorruptProcs /\ n \in CorruptDl

Rep(p, n) == [i \in 1..n |-> p]
Sched ==
  CASE Witness = "unlinkGood" ->
         \* b refetches corrupt content and reaches :145; a refetches good
         \* content up to :142; b unlinks a's file; a loads, fails, raises.
         Rep("b", 9) \o Rep("a", 8) \o <<"b", "a", "a">>
    [] Witness = "cachedIgnored" ->
         \* a caches a lifter; b overwrites with corrupt content twice and
         \* unlinks; a's second call finds no file and its download fails.
         Rep("b", 2) \o Rep("a", 8) \o Rep("b", 13) \o Rep("a", 6)
    [] OTHER -> <<>>

Init ==
  /\ pc = [p \in Procs |-> "idle"]
  /\ phase = [p \in Procs |-> 1]
  /\ failed = [p \in Procs |-> FALSE]
  /\ dl = [p \in Procs |-> 0]
  /\ calls = [p \in Procs |-> 0]
  /\ dest = [c |-> InitDest, gen |-> None]
  /\ part = [n \in Names |-> NoPart]
  /\ holds = [p \in Procs |-> FALSE]
  /\ memo = [p \in Procs |-> "none"]
  /\ res = [p \in Procs |-> "none"]
  /\ blame = [p \in Procs |-> FALSE]
  /\ clobber = FALSE
  /\ lost = FALSE
  /\ clock = 1

\* The call returns or raises: back to the caller.
Finish(p, r) ==
  /\ res' = [res EXCEPT ![p] = r]
  /\ calls' = [calls EXCEPT ![p] = @ + 1]
  /\ pc' = [pc EXCEPT ![p] = "idle"]

Start(p) ==
  /\ pc[p] = "idle" /\ calls[p] < Calls
  /\ pc' = [pc EXCEPT ![p] = "exists"]
  /\ phase' = [phase EXCEPT ![p] = 1]
  /\ res' = [res EXCEPT ![p] = "none"]
  /\ blame' = [blame EXCEPT ![p] = FALSE]
  /\ UNCHANGED <<failed, dl, calls, dest, part, holds, memo, clobber, lost>>

Exists(p) ==
  /\ pc[p] = "exists"
  /\ pc' = [pc EXCEPT ![p] = IF dest.c = "absent" THEN "mk" ELSE "load"]
  /\ UNCHANGED <<phase, failed, dl, calls, dest, part, holds, memo, res,
                 blame, clobber, lost>>

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
  /\ UNCHANGED <<phase, failed, dest, memo, clobber, lost>>

\* open(name, "wb") truncates the file, or creates it if the name is gone.
Open(p) ==
  /\ pc[p] = "open"
  /\ part' = [part EXCEPT ![Name(p)] = [c |-> "partial", own |-> p]]
  /\ pc' = [pc EXCEPT ![p] = "stream"]
  /\ UNCHANGED <<phase, failed, dl, calls, dest, holds, memo, res, blame,
                 clobber, lost>>

\* Writes go through the open descriptor: they cannot recreate a name that
\* was unlinked meanwhile.
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
            /\ pc' = [pc EXCEPT ![p] = "replace"]
            /\ UNCHANGED failed
  /\ UNCHANGED <<phase, dl, calls, dest, holds, memo, res, clobber, lost>>

Replace(p) ==
  /\ pc[p] = "replace"
  /\ pc' = [pc EXCEPT ![p] = "cleanup"]
  /\ IF Fails(p, dl[p], "replace")
       THEN /\ failed' = [failed EXCEPT ![p] = TRUE]
            /\ blame' = [blame EXCEPT ![p] = TRUE]
            /\ UNCHANGED <<dest, part, holds, lost>>
       ELSE IF part[Name(p)].c = "none"
         THEN \* FileNotFoundError from os.replace: the sibling vanished.
              /\ failed' = [failed EXCEPT ![p] = TRUE]
              /\ lost' = TRUE
              /\ UNCHANGED <<dest, part, holds, blame>>
         ELSE /\ dest' = [c |-> part[Name(p)].c, gen |-> p]
              /\ part' = [part EXCEPT ![Name(p)] = NoPart]
              /\ holds' = [holds EXCEPT ![p] = FALSE]
              /\ UNCHANGED <<failed, blame, lost>>
  /\ UNCHANGED <<phase, dl, calls, memo, res, clobber>>

Cleanup(p) ==
  /\ pc[p] = "cleanup"
  /\ part' = [part EXCEPT ![Name(p)] = NoPart]
  /\ holds' = [holds EXCEPT ![p] = FALSE]
  /\ IF failed[p]
       THEN /\ failed' = [failed EXCEPT ![p] = FALSE]
            /\ Finish(p, "dde")
       ELSE /\ pc' = [pc EXCEPT ![p] = "load"]
            /\ UNCHANGED <<failed, res, calls>>
  /\ UNCHANGED <<phase, dl, dest, memo, blame, clobber, lost>>

Load(p) ==
  /\ pc[p] = "load"
  /\ IF memo[p] # "none"
       THEN Finish(p, "ok") /\ UNCHANGED <<memo, phase>>
       ELSE IF dest.c \in {"good", "partial"}
         THEN /\ memo' = [memo EXCEPT ![p] = dest.c]
              /\ Finish(p, "ok")
              /\ UNCHANGED phase
         ELSE \* ValidationError: refetch once (:139-140), then give up (:143).
              /\ pc' = [pc EXCEPT ![p] = IF phase[p] = 1 THEN "mk" ELSE "unlink"]
              /\ phase' = [phase EXCEPT ![p] = 2]
              /\ UNCHANGED <<memo, res, calls>>
  /\ UNCHANGED <<failed, dl, dest, part, holds, blame, clobber, lost>>

Unlink(p) ==
  /\ pc[p] = "unlink"
  /\ clobber' = (clobber \/ dest.c = "good")
  /\ dest' = Absent
  /\ Finish(p, "dde")
  /\ UNCHANGED <<phase, failed, dl, part, holds, memo, blame, lost>>

Step(p) == \/ Start(p) \/ Exists(p) \/ Mk(p) \/ Open(p) \/ Stream(p)
           \/ Replace(p) \/ Cleanup(p) \/ Load(p) \/ Unlink(p)

PStep(p) ==
  /\ IF clock <= Len(Sched) THEN Sched[clock] = p ELSE TRUE
  /\ Step(p)
  /\ clock' = IF clock <= Len(Sched) THEN clock + 1 ELSE clock

AllDone == \A p \in Procs : pc[p] = "idle" /\ calls[p] = Calls

\* Stutter at the end so that termination is not reported as a deadlock.
Terminated == AllDone /\ UNCHANGED vars

Next == (\E p \in Procs : PStep(p)) \/ Terminated

\* Weak fairness on every process step. No strong fairness is assumed.
Spec == Init /\ [][Next]_vars /\ \A p \in Procs : WF_vars(PStep(p))

--------------------------------------------------------------------------
TypeOK ==
  /\ pc \in [Procs -> {"idle", "exists", "mk", "open", "stream", "replace",
                       "cleanup", "load", "unlink"}]
  /\ phase \in [Procs -> 1..2]
  /\ failed \in [Procs -> BOOLEAN]
  /\ dl \in [Procs -> 0..MaxDl]
  /\ calls \in [Procs -> 0..Calls]
  /\ dest \in [c : {"absent", "good", "corrupt", "partial"}, gen : Procs \cup {None}]
  /\ part \in [Names -> [c : {"none", "empty", "partial", "good", "corrupt"},
                         own : Procs \cup {None}]]
  /\ holds \in [Procs -> BOOLEAN]
  /\ memo \in [Procs -> {"none", "good", "partial"}]
  /\ res \in [Procs -> {"none", "ok", "dde"}]
  /\ blame \in [Procs -> BOOLEAN]
  /\ clobber \in BOOLEAN /\ lost \in BOOLEAN
  /\ clock \in 1..(Len(Sched) + 1)

\* Every sibling file belongs to exactly one live staged_path frame, and
\* every such frame still has its file.
PartOwnership ==
  /\ \A n \in Names : part[n].c # "none" =>
        /\ part[n].own \in Procs
        /\ Name(part[n].own) = n
        /\ holds[part[n].own]
  /\ \A p \in Procs : holds[p] =>
        /\ pc[p] \in {"open", "stream", "replace", "cleanup"}
        /\ part[Name(p)].c # "none"
        /\ part[Name(p)].own = p

\* os.replace always finds the sibling it was given.
PartNeverVanishes == ~lost

NoLeak == AllDone => \A n \in Names : part[n] = NoPart

\* The destination never holds a half-written file.
DestNeverPartial == dest.c # "partial"

\* chain_lifter never hands out a lifter built from anything but a complete,
\* parseable chain.
LifterGood == \A p \in Procs : memo[p] \in {"none", "good"}

\* :145 never removes a good chain.
NoGoodUnlink == ~clobber

\* A call raises only if one of its own downloads failed or was corrupt.
NoSpuriousFailure == \A p \in Procs : res[p] = "dde" => blame[p]

\* A process that already holds a lifter does not raise.
CachedLifterServed == \A p \in Procs : res[p] = "dde" => memo[p] = "none"

\* At quiescence the destination is absent or good, unless the corrupt file
\* was there before any process ran.
QuiescentClean ==
  AllDone => (dest.c \in {"absent", "good"} \/ dest.gen = None)

Termination == <>AllDone

Goal ==
  CASE Witness = "unlinkGood" ->
         clobber /\ \E p \in Procs : res[p] = "dde" /\ ~blame[p]
    [] Witness = "corruptLeft" ->
         AllDone /\ dest.c = "corrupt" /\ dest.gen \in Procs
    [] Witness = "cachedIgnored" ->
         \E p \in Procs : res[p] = "dde" /\ memo[p] = "good"
    [] OTHER -> FALSE
WitnessReached == <>Goal
====
