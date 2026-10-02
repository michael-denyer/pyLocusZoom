----------------------- MODULE RecombPublishLegacySymlink -----------------------
(***************************************************************************)
(* The expected failures of RecombPublish, kept as passing witnesses.       *)
(*                                                                         *)
(* When the cache path starts as the symlink older releases published      *)
(* behind, a force writer unlinks it (recombination.py l.190) and renames  *)
(* its staging directory into place (l.196) in two steps. Between them the *)
(* path is absent, and three strict claims of RecombPublish break:         *)
(*   "gap"       NoGap: a reader that saw a complete cache (l.546) finds   *)
(*               no map at load_recombination_map's exists() (l.370).      *)
(*   "crash"     NoUndocumentedError: _holds_only_maps saw a directory     *)
(*               (l.159) and iterdir (l.160) raises FileNotFoundError.     *)
(*   "spurious"  ValidationJustified: _holds_only_maps saw the path exist  *)
(*               (l.157), is_dir (l.159) is now False, and the caller      *)
(*               raises ValidationError (l.336) for a path that held only  *)
(*               maps.                                                     *)
(*                                                                         *)
(* An invariant cannot state that a bad state is reachable, so each        *)
(* failure is pinned as a schedule instead: step i is a RecombPublish      *)
(* process step of writer Schedule[i], from one fixed initial state. Every *)
(* behaviour of this spec is therefore a behaviour of RecombPublish, and   *)
(* Witnessed passes exactly while the failure is still reachable. A change *)
(* that closes the window makes the run fail, which is the signal to move  *)
(* the claim back into the strict runs of RecombPublish.matrix.            *)
(***************************************************************************)
EXTENDS Naturals, Sequences

CONSTANTS Writers, Readers, Names, Stray, Other, InitKinds, Modes, FailAt,
          ShortAt, Sticky, ExtSteps,
          W1, W2,    \* the two writers of the schedule
          Witness    \* "gap", "crash" or "spurious"

VARIABLES out, mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname, swapped,
          foreign0, ext,
          i          \* next schedule position

M == INSTANCE RecombPublish

wvars == <<out, mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname, swapped,
           foreign0, ext, i>>

\* W1 is the force writer that replaces the symlink.
W1Mode == IF Witness = "gap" THEN "ensure" ELSE "force"
W2Mode == IF Witness = "gap" THEN "force" ELSE "forcechecked"

Schedule ==
    CASE Witness = "gap" ->
           \* W2: l.333 force. W1: l.546 exists, glob, cache hit.
           \* W2: l.285 mkdtemp, l.288 download, l.293 stage, l.180, l.187,
           \* l.190 unlink. W1: l.370 exists() finds nothing.
           <<W2, W1, W1, W1, W2, W2, W2, W2, W2, W2, W1>>
      [] Witness = "crash" ->
           \* W1: l.333, l.285, l.288, l.293. W2: l.333, l.157 exists,
           \* l.159 is_dir. W1: l.180, l.187, l.190 unlink. W2: l.160 iterdir.
           <<W1, W1, W1, W1, W2, W2, W2, W1, W1, W1, W2>>
      [] Witness = "spurious" ->
           \* W1: l.333, l.285, l.288, l.293. W2: l.333, l.157 exists.
           \* W1: l.180, l.187, l.190 unlink. W2: l.159 is_dir is False.
           <<W1, W1, W1, W1, W2, W2, W1, W1, W1, W2>>

Reached ==
    CASE Witness = "gap" -> pc[W1] = "notfound"
      [] Witness = "crash" -> pc[W2] = "crashed"
      [] Witness = "spurious" -> pc[W2] = "invalid" /\ ~foreign0

Init ==
    /\ M!Init
    /\ out = [kind |-> "symlink", files |-> Names]
    /\ mode = [w \in Writers |-> IF w = W1 THEN W1Mode ELSE W2Mode]
    /\ rname = [w \in Writers |-> M!Min(Names)]
    /\ i = 1

Next ==
    \/ /\ i <= Len(Schedule)
       /\ M!Proc(Schedule[i])
       /\ i' = i + 1
    \/ /\ i > Len(Schedule)
       /\ UNCHANGED wvars

Spec == Init /\ [][Next]_wvars /\ WF_wvars(Next)

TypeOK == M!TypeOK /\ i \in 1..(Len(Schedule) + 1)

\* The schedule runs to its end and ends in the failure.
Witnessed == <>(i > Len(Schedule) /\ Reached)
=============================================================================
