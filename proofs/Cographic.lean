import Std

namespace IndependentSAT

def equalEdge (a b : Bool) : Nat :=
  if a = b then 1 else 0

theorem split_triangle_iff (a b c : Bool) :
    equalEdge a b + equalEdge a c + equalEdge b c = 1 ↔
      ¬ (a = b ∧ b = c) := by
  cases a <;> cases b <;> cases c <;> decide

theorem apex_triangle_iff (a b : Bool) :
    equalEdge true a + equalEdge true b + equalEdge a b = 1 ↔
      ¬ (a = true ∧ b = true) := by
  cases a <;> cases b <;> decide

def Diamond (a b p q : Bool) : Prop :=
  a.toNat + p.toNat + q.toNat = 1 ∧
    b.toNat + p.toNat + q.toNat = 1

theorem diamond_iff (a b : Bool) :
    (∃ p q : Bool, Diamond a b p q) ↔ a = b := by
  constructor
  · rintro ⟨p, q, h⟩
    obtain ⟨ha, hb⟩ := h
    have hbits : a.toNat = b.toNat := by omega
    cases a <;> cases b <;> simp_all
  · intro h
    subst b
    cases a
    · exact ⟨true, false, by decide⟩
    · exact ⟨false, false, by decide⟩

#print axioms split_triangle_iff
#print axioms apex_triangle_iff
#print axioms diamond_iff

end IndependentSAT
