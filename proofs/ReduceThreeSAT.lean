import LocalGadget

namespace IndependentSAT.Reduction

abbrev Triple (α : Type) := α × α × α
abbrev Literal (α : Type) := α × Bool
abbrev Clause (α : Type) := Triple (Literal α)

inductive Slot where
  | first | second | third

inductive Constant where
  | zero | zeroTwin | pairLeft | pairRight | one

/-- A clause itself names its auxiliaries, so repeated clauses can share a gadget. -/
inductive Variable (α : Type) where
  | literal : α → Bool → Variable α
  | constant : Constant → Variable α
  | selector : Clause α → Slot → Variable α
  | slack : Clause α → Slot → Variable α

def literalValue (assignment : α → Bool) (lit : Literal α) : Bool :=
  if lit.2 then assignment lit.1 else !(assignment lit.1)

def clauseHolds (assignment : α → Bool) (c : Clause α) : Prop :=
  literalValue assignment c.1 = true ∨
    literalValue assignment c.2.1 = true ∨ literalValue assignment c.2.2 = true

def exact (assignment : Variable α → Bool) (row : Triple (Variable α)) : Prop :=
  (assignment row.1).toNat + (assignment row.2.1).toNat + (assignment row.2.2).toNat = 1

def exactAll (assignment : Variable α → Bool) (rows : List (Triple (Variable α))) : Prop :=
  ∀ row ∈ rows, exact assignment row

private def zeroRows : List (Triple (Variable α)) :=
  [(.constant .zero, .constant .pairLeft, .constant .pairRight),
   (.constant .zeroTwin, .constant .pairLeft, .constant .pairRight),
   (.constant .zero, .constant .zeroTwin, .constant .one)]

private def complementRow (v : α) : Triple (Variable α) :=
  (.literal v true, .literal v false, .constant .zero)

private def opposite (literal : Literal α) : Variable α :=
  .literal literal.1 (!literal.2)

private def orRows (c : Clause α) : List (Triple (Variable α)) :=
  [(.selector c .first, .selector c .second, .selector c .third),
   (.selector c .first, opposite c.1, .slack c .first),
   (.selector c .second, opposite c.2.1, .slack c .second),
   (.selector c .third, opposite c.2.2, .slack c .third)]

private def block (c : Clause α) : List (Triple (Variable α)) :=
  [complementRow c.1.1, complementRow c.2.1.1, complementRow c.2.2.1] ++ orRows c

def compile (formula : List (Clause α)) : List (Triple (Variable α)) :=
  zeroRows ++ formula.flatMap block

private theorem exact_append (g : Variable α → Bool) (xs ys : List (Triple (Variable α))) :
    exactAll g (xs ++ ys) ↔ exactAll g xs ∧ exactAll g ys := by
  simp only [exactAll, List.mem_append]
  exact ⟨fun h => ⟨fun r hr => h r (.inl hr), fun r hr => h r (.inr hr)⟩,
    fun ⟨hx, hy⟩ r hr => hr.elim (hx r) (hy r)⟩

private theorem exact_flatMap (g : Variable α → Bool) (formula : List (Clause α)) :
    exactAll g (formula.flatMap block) ↔ ∀ c ∈ formula, exactAll g (block c) := by
  simp only [exactAll, List.mem_flatMap]
  exact ⟨fun h c hc r hr => h r ⟨c, hc, hr⟩,
    fun h r ⟨c, hc, hr⟩ => h c hc r hr⟩

private theorem zero_from_rows (g : Variable α → Bool) (h : exactAll g zeroRows) :
    g (.constant .zero) = false := by
  simp only [exactAll, zeroRows, List.forall_mem_cons, List.forall_mem_nil, and_true] at h
  obtain ⟨h₁, h₂, h₃⟩ := h
  simp only [exact] at h₁ h₂ h₃
  cases hz : g (.constant .zero)
  · rfl
  · simp [hz] at h₁ h₃
    omega

private theorem opposite_of_complement (g : Variable α → Bool) (lit : Literal α)
    (hz : g (.constant .zero) = false) (h : exact g (complementRow lit.1)) :
    g (opposite lit) = !(literalValue (fun v => g (.literal v true)) lit) := by
  rcases lit with ⟨v, sign⟩
  cases sign <;> cases hp : g (.literal v true) <;> cases hn : g (.literal v false) <;>
    simp_all [exact, complementRow, opposite, literalValue, Bool.toNat]

private theorem or_rows_iff (g : Variable α → Bool) (c : Clause α) :
    exactAll g (orRows c) ↔
      OrGadget (g (opposite c.1)).toNat (g (opposite c.2.1)).toNat
        (g (opposite c.2.2)).toNat
        (g (.selector c .first)).toNat (g (.selector c .second)).toNat
        (g (.selector c .third)).toNat
        (g (.slack c .first)).toNat (g (.slack c .second)).toNat (g (.slack c .third)).toNat := by
  simp only [exactAll, orRows, List.forall_mem_cons, List.forall_mem_nil, and_true, exact, OrGadget]

private def choose (x y : Bool) : Slot → Bool
  | .first => x
  | .second => !x && y
  | .third => !x && !y

private def slotValue (a : α → Bool) (c : Clause α) : Slot → Bool
  | .first => literalValue a c.1
  | .second => literalValue a c.2.1
  | .third => literalValue a c.2.2

private def liftAssignment (a : α → Bool) : Variable α → Bool
  | .literal v sign => literalValue a (v, sign)
  | .constant .zero | .constant .zeroTwin | .constant .pairLeft => false
  | .constant .pairRight | .constant .one => true
  | .selector c s => choose (slotValue a c .first) (slotValue a c .second) s
  | .slack c s => slotValue a c s &&
      !(choose (slotValue a c .first) (slotValue a c .second) s)

private theorem choose_valid (x y z : Bool) (h : x = true ∨ y = true ∨ z = true) :
    OrGadget (!x).toNat (!y).toNat (!z).toNat
      (choose x y .first).toNat (choose x y .second).toNat (choose x y .third).toNat
      (x && !(choose x y .first)).toNat (y && !(choose x y .second)).toNat
      (z && !(choose x y .third)).toNat := by
  unfold OrGadget
  revert h
  cases x <;> cases y <;> cases z <;> decide

private theorem lift_opposite (a : α → Bool) (lit : Literal α) :
    liftAssignment a (opposite lit) = !(literalValue a lit) := by
  rcases lit with ⟨v, sign⟩
  cases sign <;> cases h : a v <;> simp [h, opposite, liftAssignment, literalValue]

private theorem lift_complement (a : α → Bool) (v : α) :
    exact (liftAssignment a) (complementRow v) := by
  cases h : a v <;> simp [h, exact, complementRow, liftAssignment, literalValue, Bool.toNat]

private theorem lift_zero (a : α → Bool) : exactAll (liftAssignment a) zeroRows := by
  simp [exactAll, zeroRows, exact, liftAssignment, Bool.toNat]

private theorem lift_block (a : α → Bool) (c : Clause α) (h : clauseHolds a c) :
    exactAll (liftAssignment a) (block c) := by
  rw [block, exact_append]
  refine ⟨?_, ?_⟩
  · simp [exactAll, lift_complement]
  · apply (or_rows_iff _ _).2
    simpa only [lift_opposite, liftAssignment, slotValue] using
      choose_valid (literalValue a c.1) (literalValue a c.2.1) (literalValue a c.2.2) h

private theorem bool_nat_le_one (x : Bool) : x.toNat ≤ 1 := by
  cases x <;> decide

theorem compile_correct (formula : List (Clause α)) :
    (∃ a : α → Bool, ∀ c ∈ formula, clauseHolds a c) ↔
      ∃ g : Variable α → Bool, exactAll g (compile formula) := by
  constructor
  · intro ⟨a, ha⟩
    refine ⟨liftAssignment a, (exact_append _ _ _).2 ⟨lift_zero a, ?_⟩⟩
    exact (exact_flatMap _ _).2 (fun c hc => lift_block a c (ha c hc))
  · intro ⟨g, hg⟩
    obtain ⟨hzero, hflat⟩ := (exact_append g _ _).1 hg
    have hz := zero_from_rows g hzero
    refine ⟨fun v => g (.literal v true), fun c hc => ?_⟩
    have hblock := (exact_flatMap g formula).1 hflat c hc
    have ⟨hcomp, hors⟩ := (exact_append g _ _).1 hblock
    have hfirst := opposite_of_complement g c.1 hz (hcomp _ (by simp))
    have hsecond := opposite_of_complement g c.2.1 hz (hcomp _ (by simp))
    have hthird := opposite_of_complement g c.2.2 hz (hcomp _ (by simp))
    have h := (or_rows_iff g c).1 hors
    have ht := (or_gadget_bool _ _ _ (bool_nat_le_one _) (bool_nat_le_one _)
      (bool_nat_le_one _)).1 ⟨_, _, _, _, _, _, h⟩
    rw [hfirst, hsecond, hthird] at ht
    simpa [clauseHolds, Bool.toNat_eq_zero] using ht

theorem compile_length (formula : List (Clause α)) :
    (compile formula).length = 3 + 7 * formula.length := by
  unfold compile
  simp only [List.length_append, zeroRows, List.length_cons, List.length_nil]
  have h : (formula.flatMap block).length = 7 * formula.length := by
    induction formula with
    | nil => rfl
    | cons c cs ih =>
      simp only [List.flatMap_cons, List.length_append, block, orRows,
        List.length_cons, List.length_nil, ih]
      omega
  omega

def distinct (row : Triple (Variable α)) : Prop :=
  row.1 ≠ row.2.1 ∧ row.1 ≠ row.2.2 ∧ row.2.1 ≠ row.2.2

theorem compile_distinct (formula : List (Clause α)) :
    ∀ row ∈ compile formula, distinct row := by
  intro row hm
  simp only [compile, List.mem_append, List.mem_flatMap] at hm
  rcases hm with h | ⟨c, _, h⟩
  · simp only [zeroRows, List.mem_cons, List.not_mem_nil, or_false] at h
    rcases h with rfl | rfl | rfl <;> simp [distinct]
  · simp only [block, orRows, List.cons_append, List.nil_append,
      List.mem_cons, List.not_mem_nil, or_false] at h
    rcases h with rfl | rfl | rfl | rfl | rfl | rfl | rfl <;>
      simp [distinct, complementRow, opposite]

#print axioms compile_correct
#print axioms compile_length
#print axioms compile_distinct

end IndependentSAT.Reduction
