----------------------- MODULE RecombPublishLegacySymlink -----------------------
(***************************************************************************)
(* The expected failure of RecombPublish, kept as a passing witness.        *)
(*                                                                         *)
(* When the cache path starts as the symlink older releases published      *)
(* behind, a force writer unlinks it (recombination.py l.202) and renames  *)
(* its staging directory into place (l.208) in two steps. Between them the *)
(* path is absent, and one strict claim of RecombPublish breaks:           *)
(*   "gap"       NoGap: a reader that saw a complete cache (l.558) finds   *)
(*               no map at load_recombination_map's exists() (l.382).      *)
(*                                                                         *)
(* An invariant cannot state that a bad state is reachable, so the         *)
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
          W1, W2     \* the two writers of the schedule

VARIABLES out, mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname, swapped,
          foreign0, ext,
          i          \* next schedule position

M == INSTANCE RecombPublish

wvars == <<out, mode, pc, outcome, tmp, stgAt, stg, todo, seen, rname, swapped,
           foreign0, ext, i>>

\* W2 is the force writer that replaces the symlink; W1 is the reader.
\* W2: l.345 force. W1: l.558 exists, glob, cache hit.
\* W2: l.297 mkdtemp, l.300 download, l.305 stage, l.192, l.199,
\* l.202 unlink. W1: l.382 exists() finds nothing.
Schedule == <<W2, W1, W1, W1, W2, W2, W2, W2, W2, W2, W1>>

Reached == pc[W1] = "notfound"

Init ==
    /\ M!Init
    /\ out = [kind |-> "symlink", files |-> Names]
    /\ mode = [w \in Writers |-> IF w = W1 THEN "ensure" ELSE "force"]
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
