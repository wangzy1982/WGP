﻿/*
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
		class Real
	>
	struct SolverHeapItem {
		EquationsVariable Variable;
		int PrevSplitIndex;
		Real Size;
	};

	template<
		class EquationSystem,
		class SolverHeapItem
	>
	class DefaultSolverSpliter {
	public:
		static int GetSplitIndex(EquationSystem* equation_system, SolverHeapItem* item) {
			bool b = false;
			int next_split_index = 0;
			double ves[2] = { 1E-12, 0 };
			for (int k = 0; k < 2; ++k) {
				next_split_index = item->PrevSplitIndex + 1;
				do {
					if (next_split_index == equation_system->GetVariableCount()) {
						next_split_index = 0;
					}
					if (item->Variable.Get(next_split_index).Length() > ves[k]) {
						b = true;
						break;
					}
					++next_split_index;
				} while (next_split_index != item->PrevSplitIndex + 1);
				if (b) {
					break;
				}
			}
			if (b) {
				return next_split_index;
			}
			return 0;
		}
	};

	template<
		class EquationSystem,
		class SolverHeapItem
	>
	class DefaultSolverPriority {
	public:
		static int Compare(EquationSystem* equation_system, SolverHeapItem* item1, SolverHeapItem* item2) {
			if (item1->Size < item2->Size) {
				return -1;
			}
			if (item1->Size > item2->Size) {
				return 1;
			}
			return 0;
		}
	};

	template<
		class EquationSystem,
		class EquationsVariable,
		class IntervalVector,
		class IntervalMatrix,
		class Matrix,
		class RealInterval = Interval,
		class Real = double,
		class Spliter = DefaultSolverSpliter<EquationSystem, SolverHeapItem<EquationsVariable, Real>>,
		class Priority = DefaultSolverPriority<EquationSystem, SolverHeapItem<EquationsVariable, Real>>
	>
	class Solver {
	public:
		Solver() {
			m_max_fuzzy_root_count = 128;
			m_root_is_dirty = true;
			m_equation_system = nullptr;
		}

		void SetMaxFuzzyRootCount(int max_fuzzy_root_count) {
			m_max_fuzzy_root_count = max_fuzzy_root_count;
			m_clear_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetEquationSystem(EquationSystem* equation_system) {
			m_equation_system = equation_system;
			m_clear_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetInitialVariables(const Array<EquationsVariable>& initial_variables) {
			m_initial_variables = initial_variables;
			m_clear_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		void SetInitialVariable(const EquationsVariable& initial_variable) {
			m_initial_variables.Clear();
			m_initial_variables.Append(initial_variable);
			m_clear_roots.Clear();
			m_fuzzy_roots.Clear();
			m_root_is_dirty = true;
		}

		bool CheckInitialVariable() {
			if (!m_equation_system) {
				return true;
			}
			for (int i = 0; i < m_initial_variables.GetCount(); ++i) {
				IntervalVector value(m_equation_system->GetEquationCount());
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
				if (heap.GetCount() < m_max_fuzzy_root_count) {
					while (heap.GetCount() > 0) {
						int j = Spliter::GetSplitIndex(m_equation_system, heap.GetPointer(0));
						EquationsVariable variable1;
						EquationsVariable variable2;
						heap.GetPointer(0)->Variable.Split(j, variable1, variable2);
						HeapPop(heap);
						Iterate(&runtime, variable1, j, heap);
						Iterate(&runtime, variable2, j, heap);
						if (heap.GetCount() >= m_max_fuzzy_root_count) {
							break;
						}
					}
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
				if (Priority::Compare(m_equation_system, heap.GetPointer(j), heap.GetPointer(i)) != -1) {
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
				if (k < heap.GetCount() && Priority::Compare(m_equation_system, heap.GetPointer(k), heap.GetPointer(j)) == 1) {
					j = k;
				}
				if (Priority::Compare(m_equation_system, heap.GetPointer(j), heap.GetPointer(i)) != 1) {
					break;
				}
				SolverHeapItem<EquationsVariable, Real> t = heap.Get(i);
				heap.Set(i, heap.Get(j));
				heap.Set(j, t);
				i = j;
			}
		}

		enum class IteratedResult {
			Fuzzy,
			OnRoot,
			NoRoot
		};

		struct IteratorRuntime {
			IntervalVector q_value;
			IntervalMatrix df;
			Matrix m1;
			Matrix m2;
			IntervalMatrix bm;
			IntervalVector b_min_value;
			IntervalVector b_max_value;
			IntervalVector t_min_value;
			IntervalVector t_max_value;
			EquationsVariable min_variable;
			EquationsVariable max_variable;
			IntervalVector q_min_value;
			IntervalVector q_max_value;
			IntervalVector dv;
			EquationsVariable temp_variable;
		public:
			IteratorRuntime(int equation_count, int variable_count) : IteratorRuntime(equation_count, variable_count,
				equation_count < variable_count ? equation_count : variable_count) {
			}
			IteratorRuntime(int equation_count, int variable_count, int degree) : 
				q_value(equation_count),
				df(equation_count, variable_count),
				m1(degree, degree),
				m2(degree, degree),
				bm(degree, variable_count),
				b_min_value(degree),
				b_max_value(degree),
				t_min_value(equation_count),
				t_max_value(equation_count),
				min_variable(variable_count),
				max_variable(variable_count),
				q_min_value(equation_count),
				q_max_value(equation_count),
				dv(variable_count),
				temp_variable(variable_count) {
			}
		};

		IteratedResult Iterate(IteratorRuntime* runtime, EquationsVariable* variable, Real& size) {
			Real prev_size = -10000;
			bool slow = true;
			while (true) {
				m_equation_system->CalculateValue(*variable, runtime->q_value);
				size = 0;
				Real accept_size = 0;
				for (int i = 0; i < m_equation_system->GetEquationCount(); ++i) {
					RealInterval v = *runtime->q_value.Get(i);
					Real value_epsilon = m_equation_system->GetValueEpsilon(i);
					if (!v.IsIntersected(0, value_epsilon)) {
						size = 0;
						return IteratedResult::NoRoot;
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
					return IteratedResult::OnRoot;
				}
				if (prev_size > 0) {
					Real delta_size = prev_size - size;
					if (delta_size <= 0.1) {
						return IteratedResult::Fuzzy;
					}
					if (slow && delta_size / prev_size <= 0.01) {
						return IteratedResult::Fuzzy;
					}
				}
				prev_size = size;
				m_equation_system->CalculatePartialDerivative(*variable, runtime->df);
				bool recheck_value = false;
				bool use_default_transform = false;
				m_equation_system->Transform(*variable, runtime->q_value, runtime->df, recheck_value, use_default_transform);
				if (recheck_value) {
					for (int i = 0; i < m_equation_system->GetEquationCount(); ++i) {
						RealInterval v = *runtime->q_value.Get(i);
						Real value_epsilon = m_equation_system->GetValueEpsilon(i);
						if (!v.IsIntersected(0, value_epsilon)) {
							size = 0;
							m_equation_system->Restore();
							return IteratedResult::NoRoot;
						}
					}
				}
				slow = true;
				bool stop = true;
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					runtime->min_variable = *variable;
					runtime->max_variable = *variable;
					runtime->min_variable.Set(i, Interval(runtime->min_variable.Get(i).Min, runtime->min_variable.Get(i).Min));
					runtime->max_variable.Set(i, Interval(runtime->max_variable.Get(i).Max, runtime->max_variable.Get(i).Max));
					m_equation_system->CalculateValue(runtime->min_variable, runtime->q_min_value);
					m_equation_system->CalculateValue(runtime->max_variable, runtime->q_max_value);
					Real delta_min = 0;
					Real delta_max = 0;
					for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
						RealInterval* min_value = runtime->q_min_value.Get(j);
						RealInterval* max_value = runtime->q_max_value.Get(j);
						RealInterval* dvi = runtime->df.Get(j, i);
						if (min_value->Min > m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Min >= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = -min_value->Min / dvi->Min;
							if (d > delta_min) {
								delta_min = d;
							}
						}
						else if (min_value->Max < -m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Max <= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = -min_value->Max / dvi->Max;
							if (d > delta_min) {
								delta_min = d;
							}
						}
						if (max_value->Min > m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Max <= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = max_value->Min / dvi->Max;
							if (d > delta_max) {
								delta_max = d;
							}
						}
						else if (max_value->Max < -m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Min >= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = max_value->Max / dvi->Min;
							if (d > delta_max) {
								delta_max = d;
							}
						}
					}
					if (delta_min > 0 || delta_max > 0) {
						stop = false;
						Real d1 = variable->Get(i).Min + delta_min;
						Real d2 = variable->Get(i).Max - delta_max;
						if (d1 > d2) {
							if (variable->Get(i).Min == variable->Get(i).Max) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real dt = (variable->Get(i).Max - variable->Get(i).Min) / (delta_max + delta_min);
							double d = variable->Get(i).Min + delta_min * dt;
							variable->Set(i, Interval(d, d));
							slow = false;
						}
						else {
							if (d2 - d1 <= (variable->Get(i).Max - variable->Get(i).Min) * 0.9) {
								slow = false;
							}
							variable->Set(i, Interval(d1, d2));
						}
					}
				}
				if (!slow) {
					m_equation_system->Restore();
					continue;
				}
				if (!use_default_transform) {
					m_equation_system->Restore();
					continue;
				}
				int degree = runtime->m1.GetRowCount();
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < degree; ++j) {
						*runtime->m1.Get(i, j) = runtime->df.Get(i, j)->Center();
					}
				}
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < degree; ++j) {
						*runtime->m2.Get(i, j) = 0;
					}
					*runtime->m2.Get(i, i) = 1;
				}
				runtime->m1.GaussianElimination(runtime->m2);
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < m_equation_system->GetVariableCount(); ++j) {
						*runtime->bm.Get(i, j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->bm.Get(i, j) = *runtime->bm.Get(i, j) + *runtime->m2.Get(i, k) * *runtime->df.Get(k, j);
						}
					}
				}
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					variable->Min(runtime->temp_variable);
					m_equation_system->CalculateValue(runtime->temp_variable, runtime->t_min_value);
					variable->Max(runtime->temp_variable);
					m_equation_system->CalculateValue(runtime->temp_variable, runtime->t_max_value);
					for (int j = 0; j < degree; ++j) {
						*runtime->b_min_value.Get(j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->b_min_value.Get(j) = *runtime->b_min_value.Get(j) + *runtime->m2.Get(j, k) * *runtime->t_min_value.Get(k);
						}
						*runtime->b_max_value.Get(j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->b_max_value.Get(j) = *runtime->b_max_value.Get(j) + *runtime->m2.Get(j, k) * *runtime->t_max_value.Get(k);
						}
					}
					Real delta_min = 0;
					Real delta_max = 0;
					for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
						RealInterval min_value;
						RealInterval max_value;
						if (i == j && i < degree) {
							min_value = *runtime->b_min_value.Get(j);
							max_value = *runtime->b_max_value.Get(j);
							runtime->bm.Row(j, runtime->dv);
							for (int k = 0; k < m_equation_system->GetVariableCount(); ++k) {
								if (k == i) {
									continue;
								}
								min_value = min_value + RealInterval(0, variable->Get(k).Length()) * *runtime->dv.Get(k);
								max_value = max_value + RealInterval(-variable->Get(k).Length(), 0) * *runtime->dv.Get(k);
							}
						}
						else {
							min_value = *runtime->t_min_value.Get(j);
							max_value = *runtime->t_max_value.Get(j);
							runtime->df.Row(j, runtime->dv);
							for (int k = 0; k < degree; ++k) {
								if (k == i) {
									continue;
								}
								Real c = runtime->dv.Get(k)->Center();
								if (abs(c) > g_double_epsilon) {
									min_value = min_value - *runtime->b_min_value.Get(k) * c;
									max_value = max_value - *runtime->b_max_value.Get(k) * c;
									for (int u = 0; u < m_equation_system->GetVariableCount(); ++u) {
										*runtime->dv.Get(u) = *runtime->dv.Get(u) - *runtime->bm.Get(k, u) * c;
									}
								}
							}
							for (int k = 0; k < m_equation_system->GetVariableCount(); ++k) {
								if (k == i) {
									continue;
								}
								min_value = min_value + RealInterval(0, variable->Get(k).Length()) * *runtime->dv.Get(k);
								max_value = max_value + RealInterval(-variable->Get(k).Length(), 0) * *runtime->dv.Get(k);
							}
						}
						RealInterval* dvi = runtime->dv.Get(i);
						if (min_value.Min > m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Min >= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = -min_value.Min / dvi->Min;
							if (d > delta_min) {
								delta_min = d;
							}
						}
						else if (min_value.Max < -m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Max <= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = -min_value.Max / dvi->Max;
							if (d > delta_min) {
								delta_min = d;
							}
						}
						if (max_value.Min > m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Max <= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = max_value.Min / dvi->Max;
							if (d > delta_max) {
								delta_max = d;
							}
						}
						else if (max_value.Max < -m_equation_system->GetValueEpsilon(j)) {
							if (dvi->Min >= 0) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real d = max_value.Max / dvi->Min;
							if (d > delta_max) {
								delta_max = d;
							}
						}
					}
					if (delta_min > 0 || delta_max > 0) {
						stop = false;
						Real d1 = variable->Get(i).Min + delta_min;
						Real d2 = variable->Get(i).Max - delta_max;
						if (d1 > d2) {
							if (variable->Get(i).Min == variable->Get(i).Max) {
								m_equation_system->Restore();
								return IteratedResult::NoRoot;
							}
							Real dt = (variable->Get(i).Max - variable->Get(i).Min) / (delta_max + delta_min);
							double d = variable->Get(i).Min + delta_min * dt;
							variable->Set(i, Interval(d, d));
							slow = false;
						}
						else {
							if (d2 - d1 <= (variable->Get(i).Max - variable->Get(i).Min) * 0.9) {
								slow = false;
							}
							variable->Set(i, Interval(d1, d2));
						}
					}
				}
				if (stop) {
					m_equation_system->Restore();
					return IteratedResult::Fuzzy;
				}
			}
		}
		/*
		IteratedResult Iterate(IteratorRuntime* runtime, IntervalVector* variable, Real& size) {
			Real prev_size = -10000;
			while (true) {
				m_equation_system->CalculateValue(*variable, runtime->q_value);
				variable->Center(runtime->middle_variable);
				m_equation_system->CalculateValue(runtime->middle_variable, runtime->middle_value);
				m_equation_system->CalculatePartialDerivative(*variable, runtime->df);
				size = 0;
				Real accept_size = 0;
				for (int i = 0; i < m_equation_system->GetEquationCount(); ++i) {
					RealInterval v = *runtime->middle_value.Get(i);
					for (int j = 0; j < m_equation_system->GetVariableCount(); ++j) {
						v = v + *runtime->df.Get(i, j) * (*variable->Get(j) - *runtime->middle_variable.Get(j));
					}
					v.Intersect(*runtime->q_value.Get(i));
					Real value_epsilon = m_equation_system->GetValueEpsilon(i);
					if (!v.IsIntersected(0, value_epsilon)) {
						size = 0;
						return IteratedResult::NoRoot;
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
					return IteratedResult::OnRoot;
				}
				if (prev_size > 0) {
					Real delta_size = prev_size - size;
					if (delta_size <= 0.1) {
						return IteratedResult::Fuzzy;
					}
					if (delta_size / prev_size <= 0.01) {
						return IteratedResult::Fuzzy;
					}
				}
				prev_size = size;
				int degree = runtime->m1.GetRowCount();
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < degree; ++j) {
						*runtime->m1.Get(i, j) = runtime->df.Get(i, j)->Center();
					}
				}
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < degree; ++j) {
						*runtime->m2.Get(i, j) = 0;
					}
					*runtime->m2.Get(i, i) = 1;
				}
				runtime->m1.GaussianElimination(runtime->m2);
				for (int i = 0; i < degree; ++i) {
					for (int j = 0; j < m_equation_system->GetVariableCount(); ++j) {
						*runtime->bm.Get(i, j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->bm.Get(i, j) = *runtime->bm.Get(i, j) + *runtime->m2.Get(i, k) * *runtime->df.Get(k, j);
						}
					}
				}
				bool stop = true;
				for (int i = 0; i < m_equation_system->GetVariableCount(); ++i) {
					variable->Min(runtime->temp_variable);
					m_equation_system->CalculateValue(runtime->temp_variable, runtime->t_min_value);
					variable->Max(runtime->temp_variable);
					m_equation_system->CalculateValue(runtime->temp_variable, runtime->t_max_value);
					runtime->min_variable = *variable;
					runtime->max_variable = *variable;
					runtime->min_variable.Get(i)->Max = runtime->min_variable.Get(i)->Min;
					runtime->max_variable.Get(i)->Min = runtime->max_variable.Get(i)->Max;
					m_equation_system->CalculateValue(runtime->min_variable, runtime->q_min_value);
					m_equation_system->CalculateValue(runtime->max_variable, runtime->q_max_value);
					for (int j = 0; j < degree; ++j) {
						*runtime->b_min_value.Get(j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->b_min_value.Get(j) = *runtime->b_min_value.Get(j) + *runtime->m2.Get(j, k) * *runtime->t_min_value.Get(k);
						}
						*runtime->b_max_value.Get(j) = 0;
						for (int k = 0; k < degree; ++k) {
							*runtime->b_max_value.Get(j) = *runtime->b_max_value.Get(j) + *runtime->m2.Get(j, k) * *runtime->t_max_value.Get(k);
						}
					}
					Real delta_min = 0;
					Real delta_max = 0;
					for (int j = 0; j < m_equation_system->GetEquationCount(); ++j) {
						for (int st = 0; st <= 1; ++st) {
							RealInterval min_value;
							RealInterval max_value;
							if (st == 0) {
								min_value = *runtime->t_min_value.Get(j);
								max_value = *runtime->t_max_value.Get(j);
								runtime->dv = runtime->df.Row(j);
								for (int k = 0; k < m_equation_system->GetVariableCount(); ++k) {
									if (k == i) {
										continue;
									}
									min_value = min_value + RealInterval(0, variable->Get(k)->Length()) * *runtime->dv.Get(k);
									max_value = max_value + RealInterval(-variable->Get(k)->Length(), 0) * *runtime->dv.Get(k);
								}
								min_value.Intersect(*runtime->q_min_value.Get(j));
								max_value.Intersect(*runtime->q_max_value.Get(j));
							}
							else {
								if (i == j && i < degree) {
									min_value = *runtime->b_min_value.Get(j);
									max_value = *runtime->b_max_value.Get(j);
									runtime->dv = runtime->bm.Row(j);
									for (int k = 0; k < m_equation_system->GetVariableCount(); ++k) {
										if (k == i) {
											continue;
										}
										min_value = min_value + RealInterval(0, variable->Get(k)->Length()) * *runtime->dv.Get(k);
										max_value = max_value + RealInterval(-variable->Get(k)->Length(), 0) * *runtime->dv.Get(k);
									}
								}
								else {
									min_value = *runtime->t_min_value.Get(j);
									max_value = *runtime->t_max_value.Get(j);
									runtime->dv = runtime->df.Row(j);
									for (int k = 0; k < degree; ++k) {
										if (k == i) {
											continue;
										}
										Real c = runtime->dv.Get(k)->Center();
										if (abs(c) > g_double_epsilon) {
											min_value = min_value - *runtime->b_min_value.Get(k) * c;
											max_value = max_value - *runtime->b_max_value.Get(k) * c;
											for (int u = 0; u < m_equation_system->GetVariableCount(); ++u) {
												*runtime->dv.Get(u) = *runtime->dv.Get(u) - *runtime->bm.Get(k, u) * c;
											}
										}
									}
									for (int k = 0; k < m_equation_system->GetVariableCount(); ++k) {
										if (k == i) {
											continue;
										}
										min_value = min_value + RealInterval(0, variable->Get(k)->Length()) * *runtime->dv.Get(k);
										max_value = max_value + RealInterval(-variable->Get(k)->Length(), 0) * *runtime->dv.Get(k);
									}
								}
							}
							if (min_value.Min > m_equation_system->GetValueEpsilon(j)) {
								if (runtime->dv.Get(i)->Min >= 0) {
									return IteratedResult::NoRoot;
								}
								Real d = -min_value.Min / runtime->dv.Get(i)->Min;
								if (d > delta_min) {
									delta_min = d;
								}
							}
							else if (min_value.Max < -m_equation_system->GetValueEpsilon(j)) {
								if (runtime->dv.Get(i)->Max <= 0) {
									return IteratedResult::NoRoot;
								}
								Real d = -min_value.Max / runtime->dv.Get(i)->Max;
								if (d > delta_min) {
									delta_min = d;
								}
							}
							if (max_value.Min > m_equation_system->GetValueEpsilon(j)) {
								if (runtime->dv.Get(i)->Max <= 0) {
									return IteratedResult::NoRoot;
								}
								Real d = max_value.Min / runtime->dv.Get(i)->Max;
								if (d > delta_max) {
									delta_max = d;
								}
							}
							else if (max_value.Max < -m_equation_system->GetValueEpsilon(j)) {
								if (runtime->dv.Get(i)->Min >= 0) {
									return IteratedResult::NoRoot;
								}
								Real d = max_value.Max / runtime->dv.Get(i)->Min;
								if (d > delta_max) {
									delta_max = d;
								}
							}
						}
					}
					if (delta_min > 0 || delta_max > 0) {
						stop = false;
						Real d1 = variable->Get(i)->Min + delta_min;
						Real d2 = variable->Get(i)->Max - delta_max;
						if (d1 > d2) {
							if (variable->Get(i)->Min == variable->Get(i)->Max) {
								return IteratedResult::NoRoot;
							}
							Real dt = (variable->Get(i)->Max - variable->Get(i)->Min) / (delta_max + delta_min);
							variable->Get(i)->Min += delta_min * dt;
							variable->Get(i)->Max = variable->Get(i)->Min;
						}
						else {
							variable->Get(i)->Min = d1;
							variable->Get(i)->Max = d2;
						}
					}
				}
				if (stop) {
					return IteratedResult::Fuzzy;
				}
			}
		}
		*/

		void Iterate(IteratorRuntime* runtime, EquationsVariable variable, int prev_split_index,
			Array<SolverHeapItem<EquationsVariable, Real>>& heap) {
			Real size;
			IteratedResult r = Iterate(runtime, &variable, size);
			if (r == IteratedResult::Fuzzy) {
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
			else if (r == IteratedResult::OnRoot) {
				m_clear_roots.Append(variable);
			}
		}
	private:
		int m_max_fuzzy_root_count;
		EquationSystem* m_equation_system;
		Array<EquationsVariable> m_initial_variables;
		Array<EquationsVariable> m_clear_roots;
		Array<EquationsVariable> m_fuzzy_roots;
		bool m_root_is_dirty;
	};

}

#pragma warning(pop)

#endif