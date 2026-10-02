/-!
# The HTTP retry loop of pyLocusZoom

Source: `src/pylocuszoom/_http.py`
* `_retryable`      l.32-36
* `_with_retries`   l.45-68
* `request_json`    l.90-149 (only the exception split at l.132-144)

The model checks the code as written. It keeps the source names.

Modelling choices
* The environment is a finite list of per-attempt outcomes. `attempt()` number
  `i` (1-based) sees the `i`-th outcome. If the list runs out while the loop
  still wants another attempt, the run ends in `Exit.starved`; that is an
  artefact of the finite environment, not a behaviour of the code.
* `max_retries` is an `Int`: Python accepts 0 and negative values.
* `attempt_number` is an `Int` (Python int, unbounded, compared with
  `max_retries`).
* `delay` is a Python float that starts at `retry_delay` and is doubled. The
  model counts it in `Nat` units of `retry_delay` (1, 2, 4, ...). This assumes
  `retry_delay >= 0` and that float doubling is exact (true until overflow to
  `inf` near 2^1024). A negative `retry_delay` makes `time.sleep` raise
  `ValueError`, which is outside the model.
* `sleeps` and `calls` are ghost fields: the list of `time.sleep` arguments so
  far and the number of `attempt()` calls so far.
* Only `requests.RequestException` is modelled. Any other exception leaves the
  loop at once (l.62 does not catch it).
-/

/-- A `requests.RequestException`. `http status` is a `requests.HTTPError`;
`status` is `none` when the error carries no response (`_status_of`, l.39-42).
`connection` stands for every `RequestException` that is not an `HTTPError`. -/
inductive Err where
  | connection
  | http (status : Option Nat)
deriving Repr, DecidableEq

/-- `_http.py:27`. -/
def RETRYABLE_STATUS : List Nat := [429, 503]

/-- `_http.py:32-36`. -/
def _retryable : Err → Bool
  | .connection => true
  | .http none => false
  | .http (some status) => RETRYABLE_STATUS.contains status

/-- What one call of `attempt()` does. -/
inductive Outcome where
  | success
  | error (e : Err)
deriving Repr, DecidableEq

/-- The outcome is an error that `_retryable` accepts. -/
def Outcome.retry : Outcome → Bool
  | .success => false
  | .error e => _retryable e

structure State where
  /-- `attempt_number`, l.58 and l.68. -/
  attempt_number : Int
  /-- `delay` in units of `retry_delay`, l.57 and l.67. -/
  delay : Nat
  /-- Ghost: arguments of `time.sleep` so far, oldest first (l.66). -/
  sleeps : List Nat
  /-- Ghost: number of `attempt()` calls so far (l.61). -/
  calls : Nat
deriving Repr, DecidableEq

/-- How the loop ends. `call` is the 1-based index of the `attempt()` call that
returned or whose exception is re-raised. -/
inductive Exit where
  | returned (call : Nat)
  | raised (call : Nat) (e : Err)
  | starved
deriving Repr, DecidableEq

/-- The `while True` loop, `_http.py:59-68`. One list element is one pass. -/
def loop (max_retries : Int) (s : State) : List Outcome → Exit × State
  | [] => (.starved, s)
  -- l.61 `return attempt()`
  | .success :: _ => (.returned (s.calls + 1), { s with calls := s.calls + 1 })
  | .error e :: rest =>
    -- l.63 `if attempt_number >= max_retries or not _retryable(e): raise`
    if s.attempt_number ≥ max_retries ∨ _retryable e = false then
      (.raised (s.calls + 1) e, { s with calls := s.calls + 1 })
    else
      -- l.66 `time.sleep(delay)`, l.67 `delay *= 2`, l.68 `attempt_number += 1`
      loop max_retries
        { attempt_number := s.attempt_number + 1
          delay := s.delay * 2
          sleeps := s.sleeps ++ [s.delay]
          calls := s.calls + 1 } rest

/-- `_http.py:57-58`: `delay = retry_delay` (one unit), `attempt_number = 1`. -/
def init : State := { attempt_number := 1, delay := 1, sleeps := [], calls := 0 }

/-- `_http.py:45-68`. -/
def _with_retries (max_retries : Int) (env : List Outcome) : Exit × State :=
  loop max_retries init env

/-- How `request_json` ends, `_http.py:132-147`. -/
inductive JsonExit where
  /-- l.147: the response is decoded (JSON decoding itself is not modelled). -/
  | json (call : Nat)
  /-- l.139-140: `except requests.HTTPError`; the message has no attempt count. -/
  | httpBranch (call : Nat) (status : Option Nat)
  /-- l.141-144: the message says "failed after {reported} attempts". -/
  | countBranch (reported : Int) (call : Nat)
  | starved
deriving Repr, DecidableEq

/-- `_http.py:132-144`. `reported` is `max_retries`, exactly as l.143 formats it. -/
def request_json (max_retries : Int) (env : List Outcome) : JsonExit × State :=
  match _with_retries max_retries env with
  | (.returned k, s) => (.json k, s)
  | (.raised k (.http status), s) => (.httpBranch k status, s)
  | (.raised k .connection, s) => (.countBranch max_retries k, s)
  | (.starved, s) => (.starved, s)

/-! ## Specification -/

/-- The attempt bound `max(1, max_retries)`. -/
def bound (max_retries : Int) : Nat := if max_retries ≤ 1 then 1 else max_retries.toNat

/-- The backoff schedule: `n` sleeps of 1, 2, 4, ... units. -/
def pows (n : Nat) : List Nat := (List.range n).map (2 ^ ·)

/-- An independent oracle for the bounded check. It does not walk the list: it
finds the first outcome that is not a retryable error, clips it to the bound,
and reads off the exit, the number of calls and the sleeps. It classifies
errors itself, without `_retryable`. -/
def oracle (max_retries : Int) (env : List Outcome) : Exit × Nat × List Nat :=
  let decisive : Outcome → Bool := fun o =>
    !(o == .error .connection || o == .error (.http (some 429)) || o == .error (.http (some 503)))
  let k := min (env.findIdx decisive) ((max 1 max_retries).toNat - 1)
  match env[k]? with
  | none => (.starved, env.length, pows env.length)
  | some .success => (.returned (k + 1), k + 1, pows k)
  | some (.error e) => (.raised (k + 1) e, k + 1, pows k)

/-- Properties (a)-(d) for one run: exit, call count and sleeps match the oracle. -/
def RunOK (max_retries : Int) (env : List Outcome) : Prop :=
  let r := _with_retries max_retries env
  (r.1, r.2.calls, r.2.sleeps) = oracle max_retries env

instance (m : Int) (env : List Outcome) : Decidable (RunOK m env) := by
  unfold RunOK; infer_instance

/-- Property (e) for one run: when l.143 reports an attempt count, the count is
the number of `attempt()` calls made. -/
def MessageAccurate (max_retries : Int) (env : List Outcome) : Prop :=
  match (request_json max_retries env).1 with
  | .countBranch reported call => reported = (call : Int)
  | _ => True

instance (m : Int) (env : List Outcome) : Decidable (MessageAccurate m env) := by
  unfold MessageAccurate; split <;> infer_instance

/-! ## Bounded exhaustive check -/

def alphabet : List Outcome :=
  [.success, .error .connection, .error (.http (some 429)), .error (.http (some 503)),
   .error (.http (some 404)), .error (.http none)]

def listsOfLen : Nat → List (List Outcome)
  | 0 => [[]]
  | n + 1 => (listsOfLen n).flatMap fun l => alphabet.map (· :: l)

/-- Every outcome list over `alphabet` of length at most `n`. -/
def envsUpTo (n : Nat) : List (List Outcome) := (List.range (n + 1)).flatMap listsOfLen

/-- The integers `lo, lo + 1, ..., hi`. -/
def intsFrom (lo hi : Int) : List Int :=
  (List.range (hi - lo + 1).toNat).map fun (i : Nat) => lo + (i : Int)

/-- Runs with `lo ≤ max_retries ≤ hi` and at most `n` outcomes that break
`RunOK`; each entry is `(max_retries, env, got, expected)`. -/
def badRuns (lo hi : Int) (n : Nat) :
    List (Int × List Outcome × (Exit × Nat × List Nat) × (Exit × Nat × List Nat)) :=
  (intsFrom lo hi).flatMap fun m =>
    (envsUpTo n).filterMap fun env =>
      if RunOK m env then none
      else
        let r := _with_retries m env
        some (m, env, (r.1, r.2.calls, r.2.sleeps), oracle m env)

/-- Runs that break `MessageAccurate`; each entry is `(max_retries, env, exit)`. -/
def badMessages (lo hi : Int) (n : Nat) : List (Int × List Outcome × JsonExit) :=
  (intsFrom lo hi).flatMap fun m =>
    (envsUpTo n).filterMap fun env =>
      if MessageAccurate m env then none else some (m, env, (request_json m env).1)

-- The shipped default is `max_retries = 3`. The bound covers -2..6 and every
-- outcome list of length ≤ 5 over six outcomes (9331 lists per `max_retries`).
#eval badRuns (-2) 6 5   -- []
#guard (badRuns (-2) 6 5).isEmpty

-- Property (e) holds for `max_retries ≥ 1` ...
#eval badMessages 1 6 5   -- []
#guard (badMessages 1 6 5).isEmpty

-- ... and fails for `max_retries ≤ 0`: the message reports 0 or a negative
-- count after one real attempt. This is the code as written (l.143), so the
-- guard pins the finding instead of failing the build.
#eval badMessages (-1) 0 1
-- [(-1, [error connection], countBranch (-1) 1), (0, [error connection], countBranch 0 1)]
#guard badMessages (-1) 0 1 =
  [(-1, [.error .connection], .countBranch (-1) 1), (0, [.error .connection], .countBranch 0 1)]

-- Three attempts sleep 1 + 2 = 3 units, and nothing after the third.
#guard (_with_retries 3 [.error .connection, .error .connection, .error .connection]).2.sleeps = [1, 2]

/-! ## Proofs for every `max_retries` and every outcome list -/

/-- `_retryable` accepts exactly connection errors, 429 and 503. -/
theorem _retryable_spec (e : Err) :
    _retryable e = true ↔ e = .connection ∨ e = .http (some 429) ∨ e = .http (some 503) := by
  cases e with
  | connection => simp [_retryable]
  | http status =>
    cases status with
    | none => simp [_retryable]
    | some st => simp [_retryable, RETRYABLE_STATUS]

theorem bound_pos (m : Int) : 0 < bound m := by
  unfold bound; split <;> omega

/-- `bound` is `max(1, max_retries)`. -/
theorem bound_eq_max (m : Int) : (bound m : Int) = max 1 m := by
  unfold bound; split <;> omega

theorem pows_length (n : Nat) : (pows n).length = n := by simp [pows]

/-- `1 + 2 + ... + 2^(n-1) = 2^n - 1`. -/
theorem pows_sum : ∀ n : Nat, (pows n).sum + 1 = 2 ^ n
  | 0 => by simp [pows]
  | n + 1 => by
    have ih := pows_sum n
    simp only [pows, List.range_succ, List.map_append, List.sum_append, List.map_cons,
      List.map_nil, List.sum_cons, List.sum_nil] at *
    rw [Nat.pow_succ]
    omega

/-- Every outcome in the list is a retryable error. -/
def AllRetry (l : List Outcome) : Prop := ∀ o ∈ l, o.retry = true

theorem allRetry_nil : AllRetry [] := by simp [AllRetry]

theorem allRetry_cons {o : Outcome} {l : List Outcome} :
    AllRetry (o :: l) ↔ o.retry = true ∧ AllRetry l := by simp [AllRetry]

/-- The loop state at the top of the loop after `n` failed, retried attempts. -/
def st (n : Nat) : State :=
  { attempt_number := (n : Int) + 1, delay := 2 ^ n, sleeps := pows n, calls := n }

/-- The state when the loop exits on the attempt after `n` retried ones: one
more call, no further sleep. -/
def fin (n : Nat) : State :=
  { attempt_number := (n : Int) + 1, delay := 2 ^ n, sleeps := pows n, calls := n + 1 }

theorem init_eq : init = st 0 := by simp [init, st, pows]

/-- One retry (l.66-68) takes `st n` to `st (n + 1)`. -/
theorem st_step (n : Nat) :
    ({ attempt_number := (st n).attempt_number + 1
       delay := (st n).delay * 2
       sleeps := (st n).sleeps ++ [(st n).delay]
       calls := (st n).calls + 1 } : State) = st (n + 1) := by
  simp [st, pows, List.range_succ, Nat.pow_succ]

theorem loop_error (m : Int) (n : Nat) (e : Err) (rest : List Outcome) :
    loop m (st n) (.error e :: rest) =
      if m ≤ (n : Int) + 1 ∨ _retryable e = false then (.raised (n + 1) e, fin n)
      else loop m (st (n + 1)) rest := by
  rw [← st_step]
  rfl

theorem loop_success (m : Int) (n : Nat) (rest : List Outcome) :
    loop m (st n) (.success :: rest) = (.returned (n + 1), fin n) := rfl

/-- The loop passes over a prefix of retryable errors that fits in the bound. -/
theorem loop_prefix (m : Int) (tail : List Outcome) :
    ∀ (pre : List Outcome) (n : Nat), AllRetry pre → n + pre.length < bound m →
      loop m (st n) (pre ++ tail) = loop m (st (n + pre.length)) tail := by
  intro pre
  induction pre with
  | nil => intro n _ _; simp
  | cons o pre ih =>
    intro n ha hl
    obtain ⟨ho, hp⟩ := allRetry_cons.1 ha
    simp only [List.length_cons] at hl
    cases o with
    | success => simp [Outcome.retry] at ho
    | error e =>
      simp only [Outcome.retry] at ho
      have hneg : ¬(m ≤ (n : Int) + 1 ∨ _retryable e = false) := by
        intro h
        rcases h with h | h
        · unfold bound at hl
          split at hl <;> omega
        · rw [ho] at h; cases h
      rw [List.cons_append, loop_error, if_neg hneg, ih (n + 1) hp (by omega)]
      have : n + 1 + pre.length = n + (pre.length + 1) := by omega
      rw [this]
      rfl

/-- What a run from `st n` can do. -/
def Fwd (m : Int) (n : Nat) (env : List Outcome) : Prop :=
  (loop m (st n) env = (.starved, st (n + env.length)) ∧ AllRetry env ∧
      n + env.length < bound m) ∨
  (∃ pre o rest, env = pre ++ o :: rest ∧ AllRetry pre ∧ n + pre.length < bound m ∧
    ((o = .success ∧
        loop m (st n) env = (.returned (n + pre.length + 1), fin (n + pre.length))) ∨
     (∃ e, o = .error e ∧ (_retryable e = false ∨ n + pre.length + 1 = bound m) ∧
        loop m (st n) env = (.raised (n + pre.length + 1) e, fin (n + pre.length)))))

theorem loop_fwd (m : Int) : ∀ (env : List Outcome) (n : Nat), n < bound m → Fwd m n env := by
  intro env
  induction env with
  | nil =>
    intro n hn
    exact Or.inl ⟨rfl, allRetry_nil, by simpa using hn⟩
  | cons o rest ih =>
    intro n hn
    cases o with
    | success =>
      exact Or.inr ⟨[], .success, rest, rfl, allRetry_nil, by simpa using hn,
        Or.inl ⟨rfl, rfl⟩⟩
    | error e =>
      unfold Fwd
      rw [loop_error]
      by_cases hc : m ≤ (n : Int) + 1 ∨ _retryable e = false
      · rw [if_pos hc]
        refine Or.inr ⟨[], .error e, rest, rfl, allRetry_nil, by simpa using hn,
          Or.inr ⟨e, rfl, ?_, rfl⟩⟩
        rcases hc with hc | hc
        · right
          unfold bound at hn ⊢
          simp only [List.length_nil, Nat.add_zero]
          split at hn <;> split <;> omega
        · exact Or.inl hc
      · rw [if_neg hc]
        have hm : ¬ m ≤ (n : Int) + 1 := fun h => hc (Or.inl h)
        have hr : _retryable e = true := by
          cases h : _retryable e with
          | true => rfl
          | false => exact absurd (Or.inr h) hc
        have hn' : n + 1 < bound m := by
          unfold bound; split <;> omega
        rcases ih (n + 1) hn' with ⟨h1, h2, h3⟩ | ⟨pre, o, rest', he, hp, hl, hres⟩
        · refine Or.inl ⟨?_, allRetry_cons.2 ⟨hr, h2⟩, ?_⟩
          · rw [h1]
            have : n + 1 + rest.length = n + (rest.length + 1) := by omega
            rw [this]; rfl
          · simp only [List.length_cons]; omega
        · have hidx : n + 1 + pre.length = n + (pre.length + 1) := by omega
          rw [hidx] at hl hres
          exact Or.inr ⟨.error e :: pre, o, rest', by rw [he]; rfl,
            allRetry_cons.2 ⟨hr, hp⟩, hl, hres⟩

/-- Master case split for `_with_retries`. Either the outcome list ran out
(all of it retryable errors, fewer than the bound), or the run stopped on the
first outcome `o` that follows a prefix `pre` of retryable errors: a success
returns, an error is re-raised because it is fatal or because this is attempt
number `bound max_retries`. In both stopping cases the final state is
`fin pre.length`: `pre.length + 1` calls and `pre.length` sleeps. -/
theorem with_retries_cases (m : Int) (env : List Outcome) :
    (_with_retries m env = (.starved, st env.length) ∧ AllRetry env ∧ env.length < bound m) ∨
    (∃ pre o rest, env = pre ++ o :: rest ∧ AllRetry pre ∧ pre.length < bound m ∧
      ((o = .success ∧ _with_retries m env = (.returned (pre.length + 1), fin pre.length)) ∨
       (∃ e, o = .error e ∧ (_retryable e = false ∨ pre.length + 1 = bound m) ∧
          _with_retries m env = (.raised (pre.length + 1) e, fin pre.length)))) := by
  have h := loop_fwd m env 0 (bound_pos m)
  simp only [Fwd, Nat.zero_add] at h
  unfold _with_retries
  rw [init_eq]
  exact h

theorem mid_getElem (pre : List Outcome) (o : Outcome) (rest : List Outcome) :
    (pre ++ o :: rest)[pre.length]? = some o := by simp

/-- (a) Attempt bound, and so termination: `attempt()` is called at most
`max(1, max_retries)` times, for every `max_retries : Int`. No hypotheses. -/
theorem calls_le_bound (m : Int) (env : List Outcome) :
    (_with_retries m env).2.calls ≤ bound m := by
  rcases with_retries_cases m env with ⟨h, _, hl⟩ | ⟨pre, o, rest, _, _, hl, ⟨_, h⟩ | ⟨e, _, _, h⟩⟩
  · rw [h]; simp only [st]; omega
  · rw [h]; simp only [fin]; omega
  · rw [h]; simp only [fin]; omega

/-- The run never asks for more outcomes than the list holds unless it is starved. -/
theorem not_starved_of_length (m : Int) (env : List Outcome) (hlen : bound m ≤ env.length) :
    (_with_retries m env).1 ≠ .starved := by
  rcases with_retries_cases m env with ⟨_, _, hl⟩ | ⟨pre, o, rest, _, _, _, ⟨_, h⟩ | ⟨e, _, _, h⟩⟩
  · omega
  · rw [h]; simp
  · rw [h]; simp

/-- (a) Exactly the bound when every outcome is a retryable error. Hypotheses:
`hall` (all outcomes retryable) and `hlen` (the environment supplies at least
`bound` outcomes). The error raised is that of attempt `bound`. -/
theorem all_retryable_exact (m : Int) (env : List Outcome) (hall : AllRetry env)
    (hlen : bound m ≤ env.length) :
    ∃ e, _with_retries m env = (.raised (bound m) e, fin (bound m - 1)) ∧
      env[bound m - 1]? = some (.error e) := by
  rcases with_retries_cases m env with
    ⟨_, _, hl⟩ | ⟨pre, o, rest, he, _, _, ⟨ho, _⟩ | ⟨e, ho, hd, h⟩⟩
  · omega
  · have hmem : o ∈ env := by rw [he]; simp
    have := hall o hmem
    rw [ho] at this
    simp [Outcome.retry] at this
  · have hmem : o ∈ env := by rw [he]; simp
    have hr := hall o hmem
    rw [ho] at hr
    simp only [Outcome.retry] at hr
    rcases hd with hd | hd
    · rw [hr] at hd; cases hd
    · refine ⟨e, ?_, ?_⟩
      · rw [h, ← hd]; rfl
      · rw [← hd, he, ← ho]
        exact mid_getElem pre o rest

/-- (b) The run returns on call `k` iff a success sits at position `k`, every
earlier outcome is a retryable error (so no fatal error comes first), and
`k ≤ bound`. No hypotheses. -/
theorem returned_iff (m : Int) (env : List Outcome) (k : Nat) :
    (_with_retries m env).1 = .returned k ↔
      ∃ pre rest, env = pre ++ .success :: rest ∧ AllRetry pre ∧ pre.length < bound m ∧
        k = pre.length + 1 := by
  constructor
  · intro hk
    rcases with_retries_cases m env with
      ⟨h, _, _⟩ | ⟨pre, o, rest, he, hp, hl, ⟨ho, h⟩ | ⟨e, _, _, h⟩⟩
    · rw [h] at hk; simp at hk
    · rw [h] at hk
      simp only [Exit.returned.injEq] at hk
      subst ho
      exact ⟨pre, rest, he, hp, hl, hk.symm⟩
    · rw [h] at hk; simp at hk
  · rintro ⟨pre, rest, rfl, hp, hl, rfl⟩
    unfold _with_retries
    rw [init_eq, loop_prefix m _ pre 0 hp (by omega), loop_success]
    simp

/-- (b) A fatal error is re-raised on the spot: after a prefix of retryable
errors that fits in the bound, the fatal error of call `pre.length + 1` is the
one raised, with `pre.length` sleeps in total, so none after it. Hypotheses:
`hp` (prefix retryable), `hf` (`_retryable e = false`), `hl` (within bound). -/
theorem fatal_immediate (m : Int) (pre rest : List Outcome) (e : Err) (hp : AllRetry pre)
    (hf : _retryable e = false) (hl : pre.length < bound m) :
    _with_retries m (pre ++ .error e :: rest) = (.raised (pre.length + 1) e, fin pre.length) := by
  unfold _with_retries
  rw [init_eq, loop_prefix m _ pre 0 hp (by omega), loop_error, if_pos (Or.inr hf)]
  simp

/-- (c) No sleep after the final attempt, and the backoff schedule: a run that
stops after `k` calls slept exactly `1, 2, ..., 2^(k-2)` units, `k - 1` sleeps.
Hypothesis: the run is not starved. -/
theorem exit_sleeps (m : Int) (env : List Outcome)
    (hne : (_with_retries m env).1 ≠ .starved) :
    (_with_retries m env).2.sleeps = pows ((_with_retries m env).2.calls - 1) ∧
      1 ≤ (_with_retries m env).2.calls := by
  rcases with_retries_cases m env with ⟨h, _, _⟩ | ⟨pre, o, rest, _, _, _, ⟨_, h⟩ | ⟨e, _, _, h⟩⟩
  · rw [h] at hne; simp at hne
  · rw [h]; simp [fin]
  · rw [h]; simp [fin]

/-- (c) Total sleep for `k` calls is `2^(k-1) - 1` units of `retry_delay`. -/
theorem total_sleep (m : Int) (env : List Outcome)
    (hne : (_with_retries m env).1 ≠ .starved) :
    (_with_retries m env).2.sleeps.sum = 2 ^ ((_with_retries m env).2.calls - 1) - 1 ∧
      (_with_retries m env).2.sleeps.length + 1 = (_with_retries m env).2.calls := by
  obtain ⟨hs, hc⟩ := exit_sleeps m env hne
  rw [hs, pows_length]
  have := pows_sum ((_with_retries m env).2.calls - 1)
  omega

/-- (c) While the loop is still running (starved), it has slept once per call. -/
theorem starved_sleeps (m : Int) (env : List Outcome)
    (h : (_with_retries m env).1 = .starved) :
    (_with_retries m env).2.sleeps = pows env.length ∧
      (_with_retries m env).2.calls = env.length := by
  rcases with_retries_cases m env with ⟨h', _, _⟩ | ⟨pre, o, rest, _, _, _, ⟨_, h'⟩ | ⟨e, _, _, h'⟩⟩
  · rw [h']; simp [st]
  · rw [h'] at h; simp at h
  · rw [h'] at h; simp at h

/-- (d) The exception raised is the last attempt's: if call `k` raises `e`,
then `k` is the total number of calls, outcome `k` of the environment is that
error, every earlier outcome was a retryable error, and the raise is due to a
fatal error or to reaching the bound. No hypotheses. -/
theorem raised_is_last (m : Int) (env : List Outcome) (k : Nat) (e : Err)
    (hk : (_with_retries m env).1 = .raised k e) :
    (_with_retries m env).2.calls = k ∧ env[k - 1]? = some (.error e) ∧
      (_retryable e = false ∨ k = bound m) ∧
      ∃ pre rest, env = pre ++ .error e :: rest ∧ AllRetry pre ∧ k = pre.length + 1 := by
  rcases with_retries_cases m env with
    ⟨h, _, _⟩ | ⟨pre, o, rest, he, hp, _, ⟨_, h⟩ | ⟨e', ho, hd, h⟩⟩
  · rw [h] at hk; simp at hk
  · rw [h] at hk; simp at hk
  · rw [h] at hk
    simp only [Exit.raised.injEq] at hk
    obtain ⟨hk1, hk2⟩ := hk
    subst hk1 hk2 ho
    refine ⟨by rw [h]; rfl, ?_, hd, pre, rest, he, hp, rfl⟩
    rw [he]
    exact mid_getElem pre _ rest

/-- (d) Likewise the value returned is the last call's. -/
theorem returned_is_last (m : Int) (env : List Outcome) (k : Nat)
    (hk : (_with_retries m env).1 = .returned k) :
    (_with_retries m env).2.calls = k ∧ env[k - 1]? = some .success := by
  rcases with_retries_cases m env with
    ⟨h, _, _⟩ | ⟨pre, o, rest, he, _, _, ⟨ho, h⟩ | ⟨e', _, _, h⟩⟩
  · rw [h] at hk; simp at hk
  · rw [h] at hk
    simp only [Exit.returned.injEq] at hk
    subst hk ho
    refine ⟨by rw [h]; rfl, ?_⟩
    rw [he]
    exact mid_getElem pre _ rest
  · rw [h] at hk; simp at hk

/-- (e) Whenever `request_json` takes the l.141-144 branch, it reports
`max_retries`, while the true number of attempts is `bound max_retries`. -/
theorem count_branch_attempts (m : Int) (env : List Outcome) (reported : Int) (call : Nat)
    (h : (request_json m env).1 = .countBranch reported call) :
    reported = m ∧ call = bound m ∧ (request_json m env).2.calls = bound m := by
  unfold request_json at h ⊢
  split at h <;> simp only [reduceCtorEq, JsonExit.countBranch.injEq] at h
  next k s heq =>
    obtain ⟨h1, h2⟩ := h
    have hk : (_with_retries m env).1 = .raised k .connection := by rw [heq]
    obtain ⟨hc, _, hd, _⟩ := raised_is_last m env k .connection hk
    have hkb : k = bound m := by
      rcases hd with hd | hd
      · simp [_retryable] at hd
      · exact hd
    rw [heq] at hc
    exact ⟨h1.symm, by omega, by simpa [hkb] using hc⟩

/-- (e) The l.143 message is accurate iff `max_retries ≥ 1`. -/
theorem message_accurate_iff (m : Int) (env : List Outcome) (reported : Int) (call : Nat)
    (h : (request_json m env).1 = .countBranch reported call) :
    reported = (call : Int) ↔ 1 ≤ m := by
  obtain ⟨h1, h2, _⟩ := count_branch_attempts m env reported call h
  subst h1 h2
  unfold bound
  split <;> omega

/-- (e) The concrete mismatch: `max_retries = 0` and one connection error
reports "0 attempts" after 1 attempt. -/
theorem message_mismatch_zero :
    (request_json 0 [.error .connection]).1 = .countBranch 0 1 := by decide

/-- (e) For `max_retries ≤ 0` every count message is wrong: it reports
`max_retries ≤ 0` after exactly one attempt. -/
theorem message_wrong_nonpos (m : Int) (hm : m ≤ 0) (env : List Outcome) (reported : Int)
    (call : Nat) (h : (request_json m env).1 = .countBranch reported call) :
    call = 1 ∧ reported ≤ 0 := by
  obtain ⟨h1, h2, _⟩ := count_branch_attempts m env reported call h
  subst h1 h2
  unfold bound
  split <;> omega

/-- (e) An `HTTPError` (429 or 503) that exhausts the retries takes the l.139
branch, whose message carries no attempt count. Hypotheses: `hhttp` (every
outcome is a retryable `HTTPError`), `hlen` (at least `bound` outcomes). -/
theorem http_exhaustion_no_count (m : Int) (env : List Outcome)
    (hhttp : ∀ o ∈ env, ∃ status, o = .error (.http status) ∧ _retryable (.http status) = true)
    (hlen : bound m ≤ env.length) :
    ∃ status, (request_json m env).1 = .httpBranch (bound m) status := by
  have hall : AllRetry env := by
    intro o ho
    obtain ⟨status, rfl, hr⟩ := hhttp o ho
    exact hr
  obtain ⟨e, h, hget⟩ := all_retryable_exact m env hall hlen
  have hmem : Outcome.error e ∈ env := List.mem_of_getElem? hget
  obtain ⟨status, heq, _⟩ := hhttp _ hmem
  cases heq
  refine ⟨status, ?_⟩
  unfold request_json
  rw [h]
