---------------------------- MODULE RecombPublish ----------------------------
(***************************************************************************)
(* Concurrent publication of the recombination map cache.                   *)
(*                                                                         *)
(* Models recombination._publish_map_generation and its callers under      *)
(* POSIX rename semantics. Each filesystem call is one atomic step; the    *)
(* checks that precede a call (is_symlink, is_dir, exists) are separate    *)
(* steps, as they are in the code.                                         *)
(*                                                                         *)
(* Writers run ensure_recomb_maps (Force = FALSE) or                       *)
(* download_canine_recombination_maps(force=True) (Force = TRUE) with a    *)
(* complete staging directory already extracted. The reader runs           *)
(* get_recombination_rate_for_region: a cache-hit check, then              *)
(* load_recombination_map's exists() check, then the open.                 *)
(***************************************************************************)
EXTENDS FiniteSets

CONSTANTS Writers, Names, Stray, InitKinds

\* Everything _has_complete_maps globs with "chr*_recomb.tsv".
Globbed == Names \cup {Stray}

Kinds == {"absent", "file", "dir", "symlink"}

VARIABLES
    out,      \* the cache path: [kind, files]; for a symlink, files are its target's
    pc,       \* writer program counter
    force,    \* writer ran with force=True
    toMove,   \* files a writer has still to os.replace in the fallback
    rpc,      \* reader program counter
    rname     \* the map the reader loads

vars == <<out, pc, force, toMove, rpc, rname>>

\* What a path lookup that follows symlinks sees.
IsDir == out.kind \in {"dir", "symlink"}
Visible == IF IsDir THEN out.files ELSE {}
Complete == IsDir /\ (out.files \cap Globbed) = Names

\* InitKinds limits the starting cache path to some of Kinds, so a run can
\* exclude the legacy symlink or a stray file.
Init ==
    /\ out \in [kind : InitKinds \cap {"absent", "file"}, files : {{}}]
               \cup [kind : InitKinds \cap {"dir", "symlink"}, files : SUBSET Globbed]
    /\ pc = [w \in Writers |-> "start"]
    /\ force \in [Writers -> BOOLEAN]
    /\ toMove = [w \in Writers |-> {}]
    /\ rpc = "start"
    /\ rname \in Names

Goto(w, label) == pc' = [pc EXCEPT ![w] = label]

\* ensure_recomb_maps / download_canine_recombination_maps: cache-hit check.
Start(w) ==
    /\ pc[w] = "start"
    /\ Goto(w, IF force[w] \/ ~Complete THEN "symcheck" ELSE "hit")
    /\ UNCHANGED <<out, force, toMove, rpc, rname>>

\* if output_path.is_symlink():
SymCheck(w) ==
    /\ pc[w] = "symcheck"
    /\ Goto(w, IF out.kind = "symlink" THEN "unlink" ELSE "rename")
    /\ UNCHANGED <<out, force, toMove, rpc, rname>>

\* output_path.unlink(): FileNotFoundError if gone, EISDIR/EPERM on a directory.
Unlink(w) ==
    /\ pc[w] = "unlink"
    /\ IF out.kind \in {"symlink", "file"}
         THEN /\ out' = [kind |-> "absent", files |-> {}]
              /\ Goto(w, "rename")
         ELSE /\ Goto(w, "recheck")
              /\ UNCHANGED out
    /\ UNCHANGED <<force, toMove, rpc, rname>>

\* except OSError: if output_path.is_symlink(): raise
Recheck(w) ==
    /\ pc[w] = "recheck"
    /\ Goto(w, IF out.kind = "symlink" THEN "failed" ELSE "rename")
    /\ UNCHANGED <<out, force, toMove, rpc, rname>>

\* os.rename(staging_dir, output_path): succeeds onto nothing or an empty
\* directory; fails onto a non-empty directory, a file, or a symlink.
Rename(w) ==
    /\ pc[w] = "rename"
    /\ IF out.kind = "absent" \/ (out.kind = "dir" /\ out.files = {})
         THEN /\ out' = [kind |-> "dir", files |-> Names]
              /\ Goto(w, "done")
         ELSE /\ Goto(w, "isdir")
              /\ UNCHANGED out
    /\ UNCHANGED <<force, toMove, rpc, rname>>

\* except OSError: if not output_path.is_dir(): raise
IsDirCheck(w) ==
    /\ pc[w] = "isdir"
    /\ IF IsDir
         THEN /\ Goto(w, "replace")
              /\ toMove' = [toMove EXCEPT ![w] = Names]
         ELSE /\ Goto(w, "failed")
              /\ UNCHANGED toMove
    /\ UNCHANGED <<out, force, rpc, rname>>

\* os.replace(staging_dir / name, output_path / name), one file per step.
Replace(w) ==
    /\ pc[w] = "replace"
    /\ IF toMove[w] = {}
         THEN /\ Goto(w, "strays")
              /\ UNCHANGED <<out, toMove>>
         ELSE \E n \in toMove[w] :
                IF IsDir
                  THEN /\ out' = [out EXCEPT !.files = @ \cup {n}]
                       /\ toMove' = [toMove EXCEPT ![w] = @ \ {n}]
                       /\ UNCHANGED pc
                  ELSE /\ Goto(w, "failed")
                       /\ UNCHANGED <<out, toMove>>
    /\ UNCHANGED <<force, rpc, rname>>

\* stray.unlink(missing_ok=True) for globbed names outside the map set.
Strays(w) ==
    /\ pc[w] = "strays"
    /\ out' = IF IsDir THEN [out EXCEPT !.files = @ \ {Stray}] ELSE out
    /\ Goto(w, "done")
    /\ UNCHANGED <<force, toMove, rpc, rname>>

\* Reader: cache hit, then map_file.exists(), then pd.read_csv opens it.
ReaderStart ==
    /\ rpc = "start"
    /\ rpc' = IF Complete THEN "exists" ELSE "skip"
    /\ UNCHANGED <<out, pc, force, toMove, rname>>

ReaderExists ==
    /\ rpc = "exists"
    /\ rpc' = IF rname \in Visible THEN "open" ELSE "missing"
    /\ UNCHANGED <<out, pc, force, toMove, rname>>

ReaderOpen ==
    /\ rpc = "open"
    /\ rpc' = IF rname \in Visible THEN "ok" ELSE "missing"
    /\ UNCHANGED <<out, pc, force, toMove, rname>>

Writer(w) ==
    Start(w) \/ SymCheck(w) \/ Unlink(w) \/ Recheck(w) \/ Rename(w)
        \/ IsDirCheck(w) \/ Replace(w) \/ Strays(w)

Terminated ==
    /\ \A w \in Writers : pc[w] \in {"hit", "done", "failed"}
    /\ rpc \in {"skip", "ok", "missing"}
    /\ UNCHANGED vars

Next ==
    \/ \E w \in Writers : Writer(w)
    \/ ReaderStart \/ ReaderExists \/ ReaderOpen
    \/ Terminated

Spec == Init /\ [][Next]_vars

-----------------------------------------------------------------------------
(* Claims from the _publish_map_generation docstring.                      *)

\* "a reader finds the old file or the new one, never a gap"
NoGap == rpc # "missing"

\* "a writer that loses the race to another still succeeds"
WriterSucceeds == \A w \in Writers : pc[w] # "failed"

\* Once every writer has returned, the cache is a complete map set that the
\* next ensure_recomb_maps accepts as a hit.
Converges ==
    (\A w \in Writers : pc[w] \in {"hit", "done"}) => Complete

\* A writer never writes through a legacy symlink into a directory the
\* package does not own ("its target is not ours").
NeverWritesThroughSymlink ==
    \A w \in Writers : pc[w] = "replace" => out.kind # "symlink"
=============================================================================
