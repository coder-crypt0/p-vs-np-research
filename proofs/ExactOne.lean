import Std

namespace IndependentSAT

/-- Each entry is the number of true occurrences in one evaluated clause. -/
def Exact (counts : List Nat) : Prop := ∀ c ∈ counts, c = 1

def OddRows (counts : List Nat) : Prop := ∀ c ∈ counts, c % 2 = 1

theorem odd_cost_lower_bound {counts : List Nat} (h : OddRows counts) :
    counts.length ≤ counts.sum := by
  induction counts with
  | nil => simp
  | cons c cs ih =>
    have hc := h c (by simp)
    have ht := ih (fun x hx => h x (by simp [hx]))
    simp only [List.length_cons, List.sum_cons]
    omega

theorem exact_iff_odd_and_min_cost (counts : List Nat) :
    Exact counts ↔ OddRows counts ∧ counts.sum = counts.length := by
  induction counts with
  | nil => simp [Exact, OddRows]
  | cons c cs ih =>
    change (∀ x ∈ c :: cs, x = 1) ↔
      (∀ x ∈ c :: cs, x % 2 = 1) ∧ (c :: cs).sum = (c :: cs).length
    simp only [List.forall_mem_cons, List.sum_cons, List.length_cons]
    constructor
    · rintro ⟨rfl, htail⟩
      obtain ⟨hodds, hcost⟩ := ih.mp htail
      exact ⟨⟨by decide, hodds⟩, by omega⟩
    · rintro ⟨⟨hc, ht⟩, hcost⟩
      have lower := odd_cost_lower_bound ht
      exact ⟨by omega, ih.mpr ⟨ht, by omega⟩⟩

/-- Evaluation keeps a shared variable consistent across all its occurrences. -/
def countsFor {α : Type} (assignment : α → Bool) (clauses : List (List α)) : List Nat :=
  clauses.map (fun row => (row.filter assignment).length)

theorem sat_iff_parity_at_budget {α : Type} (clauses : List (List α)) :
    (∃ a : α → Bool, Exact (countsFor a clauses)) ↔
    ∃ a : α → Bool, OddRows (countsFor a clauses) ∧
      (countsFor a clauses).sum = clauses.length := by
  simp only [exact_iff_odd_and_min_cost, countsFor, List.length_map]

#print axioms sat_iff_parity_at_budget

/-- The twelve-clause instance has no nonnegative integer solution.
    Each index occurs in a row summing to one, so all such entries would be Boolean. -/
theorem twelve_clauses_unsat (a b c d e f g h i j k l : Nat)
    (h0 : a + g + j = 1) (h1 : a + h + i = 1) (h2 : a + k + l = 1)
    (h3 : b + g + h = 1) (h4 : b + i + l = 1) (h5 : b + j + k = 1)
    (h6 : c + d + h = 1) (h7 : c + e + l = 1) (h8 : c + f + k = 1)
    (h9 : d + e + g = 1) (h10 : d + f + j = 1) (_h11 : e + f + i = 1) :
    False := by
  omega

#print axioms twelve_clauses_unsat

end IndependentSAT
