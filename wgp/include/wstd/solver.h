/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_SOLVER_
#define _WGP_STD_SOLVER_

#pragma warning(push)
#pragma warning(disable:26495)

#include "interval.h"
#include "array.h"

namespace wgp {

	template<
		class EquationsVariable,
		class Real = double
	>
	struct SolverHeapItem {
		EquationsVariable Variable;
		int PrevSplitIndex;
		Real Size;
	};

	enum class SolverIteratedResult {
		Fuzzy,
		OnIntervalRoot,
		OnClearRoot,
		NoRoot
	};

	template<
		class EquationSystem,
		class EquationsVariable,
		class EquationsValue,
		class EquationDv,
		class IntervalMatrix,
		class RealInterval = Interval,
		class Real = double
	>
	class Solver {
	public:
		Solver() {
			m_slow_threshold = 0.1;
			m_root_is_dirty = true;
			m_equation_system = nullptr;
		}

		void SetSlowThreshold(Real slow_threshold) {
			m_slow_threshold = slow_threshold;
		}

		void SetEquationSystem(EquationSystem* equation_system) {
			m_equation_system = equation_system;
			m_clear_roots.Clear();
			m_interval_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetInitialVariables(const Array<EquationsVariable>& initial_variables) {
			m_initial_variables = initial_variables;
			m_clear_roots.Clear();
			m_interval_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetInitialVariables(Array<EquationsVariable>&& initial_variables) {
			m_initial_variables = std::forward(initial_variables);
			m_clear_roots.Clear();
			m_interval_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetInitialVariable(const EquationsVariable& initial_variable) {
			m_initial_variables.Clear();
			m_initial_variables.Append(initial_variable);
			m_clear_roots.Clear();
			m_interval_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		bool CheckInitialVariable() {
			if (!m_equation_system) {
				return true;
			}
			for (int i = 0; i < m_initial_variables.GetCount(); ++i) {
				EquationsValue value(m_equation_system->GetEquationCount());
				m_equation_system->CalculateValue(m_initial_variables.Get(i), value);
				for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
					if (!(value.Get(j)->Length() <= 1E200)) {
						return false;
					}
				}
			}
			return true;
		}

		const Array<EquationsVariable>& GetClearRoots() {
			Execute();
			return m_clear_roots;
		}

		const Array<EquationsVariable>& GetIntervalRoots() {
			Execute();
			return m_interval_roots;
		}

		const Array<EquationsVariable>& GetFuzzyRoots() {
			Execute();
			return m_fuzzy_roots;
		}
	private:
		void Execute() {
			if (m_root_is_dirty) {
				IteratorRuntime runtime(m_equation_system->GetEquationCount(), m_equation_system->GetVariableCount());
				Array<SolverHeapItem<EquationsVariable, Real>> heap;
				for (int i = 0; i < m_initial_variables.GetCount(); ++i) {
					Iterate(&runtime, m_initial_variables.Get(i), -1, heap);
				}
				while (heap.GetCount() > 0) {
					if (m_equation_system->CheckFinished(heap)) {
						break;
					}
					SolverHeapItem<EquationsVariable, Real>* item = heap.GetPointer(0);
					int j = m_equation_system->GetSplitIndex(item->Variable, item->PrevSplitIndex, item->Size);
					EquationsVariable variable1;
					EquationsVariable variable2;
					heap.GetPointer(0)->Variable.Split(j, variable1, variable2);
					HeapPop(heap);
					Iterate(&runtime, variable1, j, heap);
					Iterate(&runtime, variable2, j, heap);
				}
				for (int i = 0; i < heap.GetCount(); ++i) {
					m_fuzzy_roots.Append(heap.GetPointer(i)->Variable);
				}
				m_root_is_dirty = false;
			}
		}

		void HeapPush(Array<SolverHeapItem<EquationsVariable, Real>>& heap, Real size, int prev_split_index, const EquationsVariable& variable) {
			SolverHeapItem<EquationsVariable, Real> item;
			item.Size = size;
			item.PrevSplitIndex = prev_split_index;
			item.Variable = variable;
			heap.Append(item);
			int i = heap.GetCount() - 1;
			while (i > 0) {
				int j = (i - 1) / 2;
				SolverHeapItem<EquationsVariable, Real>* itemi = heap.GetPointer(i);
				SolverHeapItem<EquationsVariable, Real>* itemj = heap.GetPointer(j);
				if (m_equation_system->CompareIteratePriority(itemj->Variable, itemj->Size, itemi->Variable, itemi->Size) != -1) {
					break;
				}
				SolverHeapItem<EquationsVariable, Real> t = heap.Get(i);
				heap.Set(i, heap.Get(j));
				heap.Set(j, t);
				i = j;
			}
		}

		void HeapPop(Array<SolverHeapItem<EquationsVariable, Real>>& heap) {
			if (heap.GetCount() > 1) {
				heap.Set(0, heap.Get(heap.GetCount() - 1));
			}
			heap.PopLast();
			int i = 0;
			while (true) {
				int j = i * 2 + 1;
				if (j >= heap.GetCount()) {
					break;
				}
				int k = j + 1;
				SolverHeapItem<EquationsVariable, Real>* itemk = heap.GetPointer(k);
				SolverHeapItem<EquationsVariable, Real>* itemj = heap.GetPointer(j);
				if (k < heap.GetCount() && m_equation_system->CompareIteratePriority(itemk->Variable, itemk->Size, itemj->Variable, itemj->Size) == 1) {
					j = k;
				}
				SolverHeapItem<EquationsVariable, Real>* itemi = heap.GetPointer(i);
				itemj = heap.GetPointer(j);
				if (m_equation_system->CompareIteratePriority(itemj->Variable, itemj->Size, itemi->Variable, itemi->Size) != 1) {
					break;
				}
				SolverHeapItem<EquationsVariable, Real> t = heap.Get(i);
				heap.Set(i, heap.Get(j));
				heap.Set(j, t);
				i = j;
			}
		}

		struct IteratorRuntime {
			EquationsValue q_value;
			IntervalMatrix df;
			EquationsVariable min_variable;
			EquationsVariable max_variable;
			EquationsValue q_min_value;
			EquationsValue q_max_value;
			EquationDv dv;
		public:
			IteratorRuntime(int equation_count, int variable_count) : IteratorRuntime(equation_count, variable_count,
				equation_count < variable_count ? equation_count : variable_count) {
			}
			IteratorRuntime(int equation_count, int variable_count, int degree) : 
				q_value(equation_count),
				df(equation_count, variable_count),
				min_variable(variable_count),
				max_variable(variable_count),
				q_min_value(equation_count),
				q_max_value(equation_count),
				dv(variable_count) {
			}
		};

		bool C1(double d) {
			return *(((int*)&d) + 1) & ((*(int*)&g_double_epsilon) << 31) & 0x80000000;
		}

		bool C2(double d) {
			return (*(((int*)&d) + 1) & 0x80000000) ^ ((~(*(int*)&g_double_epsilon) << 6) & 0x80000000);
		}

		bool C3(double d) {
			return *(((int*)&d) + 1) & (~(*(int*)&g_double_epsilon) << 9) & 0x80000000;
		}

		bool C4(double d) {
			return (*(((int*)&d) + 1) & 0x80000000) ^ (((*(int*)&g_double_epsilon) << 5) & 0x80000000);
		}

		bool C5(double d) {
			return *(((int*)&d) + 1) & (~(*(int*)&g_double_epsilon) << 13) & 0x80000000;
		}

		bool C6(double d) {
			return (*(((int*)&d) + 1) & 0x80000000) ^ (((*(int*)&g_double_epsilon) << 18) & 0x80000000);
		}

		SolverIteratedResult Iterate(IteratorRuntime* runtime, EquationsVariable* variable, Real& size) {
			Real prev_size = -10000;
			bool slow = false;
			while (true) {
				SolverIteratedResult result;
				if (m_equation_system->PreIterate(variable, result, size)) {
					return result;
				}
				m_equation_system->CalculateValue(*variable, runtime->q_value);
				size = 0;
				Real accept_size = 0;
				for (int i = 0; i < m_equation_system->GetEquationCount(); ++i) {
					RealInterval v = runtime->q_value.Get(i);
					Real value_epsilon = m_equation_system->GetValueEpsilon(i, true);
					if (!v.IsIntersected(0, value_epsilon)) {
						size = 0;
						return SolverIteratedResult::NoRoot;
					}
					Real d = v.Length();
					Real d1 = d / value_epsilon;
					if (d1 > size) {
						size = d1;
					}
					Real d2 = d / value_epsilon;
					if (d2 > accept_size) {
						accept_size = d2;
					}
				}
				if (accept_size <= 1) {
					return SolverIteratedResult::OnClearRoot;
				}
				if (prev_size > 0) {
					Real delta_size = prev_size - size;
					if (delta_size <= 0.1) {
						return SolverIteratedResult::Fuzzy;
					}
					if (slow && delta_size / prev_size <= m_slow_threshold) {
						return SolverIteratedResult::Fuzzy;
					}
				}
				prev_size = size;
				m_equation_system->CalculatePartialDerivative(*variable, runtime->df);
				slow = true;
				bool stop = true;
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					if (variable->Get(i).Length() <= g_double_epsilon) {
						continue;
					}
					runtime->min_variable = *variable;
					runtime->max_variable = *variable;
					runtime->min_variable.Set(i, Interval(runtime->min_variable.Get(i).Min, runtime->min_variable.Get(i).Min));
					runtime->max_variable.Set(i, Interval(runtime->max_variable.Get(i).Max, runtime->max_variable.Get(i).Max));
					m_equation_system->CalculateValue(runtime->min_variable, runtime->q_min_value);
					m_equation_system->CalculateValue(runtime->max_variable, runtime->q_max_value);
					Real delta_min = 0;
					Real delta_max = 0;
					for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
						RealInterval min_value = runtime->q_min_value.Get(j);
						RealInterval max_value = runtime->q_max_value.Get(j);
						RealInterval* dvi = runtime->df.Get(j, i);
						if (C2(min_value.Min - m_equation_system->GetValueEpsilon(j, true))) {
							if (dvi->Min >= 0) {
								return SolverIteratedResult::NoRoot;
							}
							else {
								Real d = -(min_value.Min - m_equation_system->GetValueEpsilon(j, false)) / dvi->Min;
								if (d > delta_min) {
									delta_min = d;
								}
							}
						}
						else if (C3(min_value.Max + m_equation_system->GetValueEpsilon(j, true))) {
							if (dvi->Max <= 0) {
								return SolverIteratedResult::NoRoot;
							}
							else {
								Real d = -(min_value.Max + m_equation_system->GetValueEpsilon(j, false)) / dvi->Max;
								if (d > delta_min) {
									delta_min = d;
								}
							}
						}
						if (C4(max_value.Min - m_equation_system->GetValueEpsilon(j, true))) {
							if (dvi->Max <= 0) {
								return SolverIteratedResult::NoRoot;
							}
							else {
								Real d = (max_value.Min - m_equation_system->GetValueEpsilon(j, false)) / dvi->Max;
								if (d > delta_max) {
									delta_max = d;
								}
							}
						}
						else if (C1(max_value.Max + m_equation_system->GetValueEpsilon(j, true))) {
							if (dvi->Min >= 0) {
								return SolverIteratedResult::NoRoot;
							}
							else {
								Real d = (max_value.Max + m_equation_system->GetValueEpsilon(j, false)) / dvi->Min;
								if (d > delta_max) {
									delta_max = d;
								}
							}
						}
					}
					if (delta_min > 0 || delta_max > 0) {
						stop = false;
						Real d1 = variable->Get(i).Min + delta_min;
						Real d2 = variable->Get(i).Max - delta_max;
						if (d1 > d2) {
							if (variable->Get(i).Min == variable->Get(i).Max) {
								return SolverIteratedResult::NoRoot;
							}
							Real dt = (variable->Get(i).Max - variable->Get(i).Min) / (delta_max + delta_min);
							double d = variable->Get(i).Min + delta_min * dt;
							variable->Set(i, Interval(d, d));
							slow = false;
						}
						else {
							if (d2 - d1 <= (variable->Get(i).Max - variable->Get(i).Min) * 0.5) {
								slow = false;
							}
							variable->Set(i, Interval(d1, d2));
						}
					}
				}
				if (stop) {
					return SolverIteratedResult::Fuzzy;
				}
			}
		}
		
		void Iterate(IteratorRuntime* runtime, EquationsVariable variable, int prev_split_index,
			Array<SolverHeapItem<EquationsVariable, Real>>& heap) {
			Real size;
			SolverIteratedResult r = Iterate(runtime, &variable, size);
			if (r == SolverIteratedResult::Fuzzy) {
				bool b = true;
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					if (variable.Get(i).Length() > m_equation_system->GetVariableEpsilon(i)) {
						b = false;
						break;
					}
				}
				if (b) {
					m_clear_roots.Append(variable);
				}
				else {
					HeapPush(heap, size, prev_split_index, variable);
				}
			}
			else if (r == SolverIteratedResult::OnIntervalRoot) {
				m_interval_roots.Append(variable);
			}
			else if (r == SolverIteratedResult::OnClearRoot) {
				m_clear_roots.Append(variable);
			}
		}
	private:
		Real m_slow_threshold;
		EquationSystem* m_equation_system;
		Array<EquationsVariable> m_initial_variables;
		Array<EquationsVariable> m_clear_roots;
		Array<EquationsVariable> m_interval_roots;
		Array<EquationsVariable> m_fuzzy_roots;
		bool m_root_is_dirty;
	};

	template<
		class EquationSystem,
		class EquationsVariable,
		class EquationsValue,
		class DeltaVariable,
		class Matrix,
		class RealInterval = Interval,
		class Real = double
	>
	class NewtonSolver {
	public:
		NewtonSolver() {
			m_equation_system = nullptr;
		}

		void SetEquationSystem(EquationSystem* equation_system) {
			m_equation_system = equation_system;
		}

		bool Iterate(EquationsVariable* variable, int max_iterate_count) {
			EquationsValue value(m_equation_system->GetEquationCount());
			DeltaVariable dv;
			Matrix m1;
			Matrix m2;
			for (int n = 0; n < max_iterate_count; ++n) {
				m_equation_system->CalculateValue(*variable, value);
				m_equation_system->CalculatePartialDerivative(*variable, m1);
				for (int i = 0; i < m2.GetRowCount(); ++i) {
					for (int j = 0; j < m2.GetColCount(); ++j) {
						*m2.Get(i, j) = 0;
					}
					*m2.Get(i, i) = 1;
				}
				m1.GaussianElimination(m2);
				for (int i = 0; i < m1.GetRowCount(); ++i) {
					if (*m1.Get(i, i) < 0.5) {
						return false;
					}
				}
				bool b = true;
				for (int i = 0; i < m_equation_system->GetEquationCount(); ++i) {
					if (!is_zero(value.Get(i), m_equation_system->GetValueEpsilon(i, true))) {
						b = false;
						break;
					}
				}
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					Real d = 0;
					for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
						d -= *m2.Get(i, j) * value.Get(j);
					}
					RealInterval vi = variable->GetInterval(i);
					Real d0 = variable->Get(i);
					if (d0 + d < vi.Min) {
						d = vi.Min - d0;
					} 
					else if (d0 + d > vi.Max) {
						d = vi.Max - d0;
					}
					if (b && !is_zero(d, m_equation_system->GetDeltaVariableEpsilon(i))) {
						b = false;
					}
					dv.Set(i, d);
				}
				if (b) {
					return true;
				}
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					variable->Set(i, variable->Get(i) + dv.Get(i));
				}
			}
			return false;
		}
	private:
		EquationSystem* m_equation_system;
	};

}

#pragma warning(pop)

#endif