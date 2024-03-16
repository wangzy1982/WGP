/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_PROF_
#define _WGP_STD_PROF_

#include "wbase.h"
#include <chrono>

namespace wgp {

	class Stopwatch {
	public:
		Stopwatch() {
			m_is_running = false;
			Reset();
		}

		void Reset() {
			if (m_is_running) {
				throw;
			}
			m_total_count = 0;
			m_total_time = std::chrono::high_resolution_clock::now() - std::chrono::high_resolution_clock::now();
		}

		void Resume() {
			if (m_is_running) {
				throw;
			}
			m_is_running = true;
			m_start_time = std::chrono::high_resolution_clock::now();
		}

		void Suspend() {
			if (!m_is_running) {
				throw;
			}
			++m_total_count;
			m_total_time += std::chrono::high_resolution_clock::now() - m_start_time;
			m_is_running = false;
		}

		double GetTotalTime() {
			if (m_is_running) {
				throw;
			}
			return m_total_time.count() * 1E-9;
		}

		double GetTotalCount() {
			if (m_is_running) {
				throw;
			}
			return m_total_count;
		}
	private:
		bool m_is_running;
		int m_total_count;
		std::chrono::steady_clock::time_point m_start_time;
		std::chrono::nanoseconds m_total_time;
	};
}

#endif