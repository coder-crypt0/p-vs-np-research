import Std

namespace IndependentSAT

def OrGadget {α : Type} [Add α] [One α]
    (nx ny nz p q r a b c : α) : Prop :=
  p + q + r = 1 ∧ p + nx + a = 1 ∧ q + ny + b = 1 ∧ r + nz + c = 1

theorem or_gadget_bool (nx ny nz : Nat) (hx : nx ≤ 1) (hy : ny ≤ 1) (hz : nz ≤ 1) :
    (∃ p q r a b c, OrGadget nx ny nz p q r a b c) ↔
      nx = 0 ∨ ny = 0 ∨ nz = 0 := by
  constructor
  · rintro ⟨p, q, r, a, b, c, hpqr, ha, hb, hc⟩
    omega
  · rintro (rfl | rfl | rfl)
    · exact ⟨1, 0, 0, 0, 1 - ny, 1 - nz, rfl, rfl, by omega, by omega⟩
    · exact ⟨0, 1, 0, 1 - nx, 0, 1 - nz, rfl, by omega, rfl, by omega⟩
    · exact ⟨0, 0, 1, 1 - nx, 1 - ny, 0, rfl, by omega, by omega, rfl⟩

theorem or_gadget_signed (nx ny nz : Int) :
    ∃ p q r a b c, OrGadget nx ny nz p q r a b c := by
  exact ⟨1, 0, 0, -nx, 1 - ny, 1 - nz, by decide,
    by omega, by omega, by omega⟩

/-- Six Fano lines fix the missing line to one and the four other points to zero. -/
theorem fano_six_lines_nat (a b c d e f g : Nat)
    (h1 : b + c + e = 1) (h2 : c + d + f = 1) (h3 : d + e + g = 1)
    (h4 : e + f + a = 1) (h5 : f + g + b = 1) (h6 : g + a + c = 1) :
    a = 1 ∧ b = 1 ∧ d = 1 ∧ c = 0 ∧ e = 0 ∧ f = 0 ∧ g = 0 := by
  omega

#print axioms or_gadget_bool
#print axioms or_gadget_signed
#print axioms fano_six_lines_nat

end IndependentSAT
