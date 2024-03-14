/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_ARRAY_
#define _WGP_STD_ARRAY_

#include <memory>

namespace wgp {

	template<class T>
	class Array {
	public:
		Array() {
			m_capacity = 0;
			m_items = nullptr;
			m_start = 0;
			m_count = 0;
		}

		Array(int capacity) {
			m_capacity = capacity;
			m_items = (T*)malloc(m_capacity * sizeof(T));
			m_start = 0;
			m_count = 0;
		}

		Array(const Array& right) {
			m_capacity = right.m_count;
			m_start = 0;
			m_count = right.m_count;
			m_items = (T*)malloc(m_capacity * sizeof(T));
			Copy(m_items, right.m_items + right.m_start, m_count);
		}

		Array(Array&& right) noexcept {
			m_capacity = right.m_capacity;
			m_start = right.m_start;
			m_count = right.m_count;
			m_items = right.m_items;
			right.m_capacity = 0;
			right.m_start = 0;
			right.m_count = 0;
			right.m_items = nullptr;
		}

		virtual ~Array() {
			FreeItems(m_items + m_start, m_count);
			free(m_items);
		}

		void operator=(const Array& right) {
			FreeItems(m_items + m_start, m_count);
			if (m_capacity >= right.m_count) {
				m_start = 0;
				m_count = right.m_count;
				Copy(m_items, right.m_items + right.m_start, m_count);
			}
			else {
				free(m_items);
				m_capacity = right.m_count;
				m_start = 0;
				m_count = right.m_count;
				m_items = (T*)malloc(m_capacity * sizeof(T));
				Copy(m_items, right.m_items + right.m_start, m_count);
			}
		}

		void operator=(Array&& right) {
			FreeItems(m_items + m_start, m_count);
			free(m_items);
			m_capacity = right.m_capacity;
			m_start = right.m_start;
			m_count = right.m_count;
			m_items = right.m_items;
			right.m_capacity = 0;
			right.m_start = 0;
			right.m_count = 0;
			right.m_items = nullptr;
		}

		void Clear() {
			FreeItems(m_items + m_start, m_count);
			m_start = 0;
			m_count = 0;
		}
		
		void Append(const T& item) {
			if (m_start + m_count == m_capacity) {
				if (m_start > 0) {
					Move(m_items, m_items + m_start, m_count);
					m_start = 0;
				}
				else {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, m_count);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
				}
			}
			new(m_items + (m_start + m_count)) T(item);
			m_count += 1;
		}

		void Append(T&& item) {
			if (m_start + m_count == m_capacity) {
				if (m_start > 0) {
					Move(m_items, m_items + m_start, m_count);
					m_start = 0;
				}
				else {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, m_count);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
				}
			}
			new(m_items + (m_start + m_count)) T(std::forward<T>(item));
			m_count += 1;
		}

		void Insert(int i, const T& item) {
			if (m_start == 0) {
				if (m_count + m_start == m_capacity) {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, i);
					Move(newItems + i + 1, m_items + m_start + i, m_count - i);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
					m_count += 1;
				}
				else {
					Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
					m_count += 1;
				}
			}
			else if (m_count + m_start == m_capacity || i < m_count * 0.5) {
				Move(m_items + m_start + i - 1, m_items + m_start + i, i);
				m_start -= 1;
				m_count += 1;
			}
			else {
				Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
				m_count += 1;
			}
			new(m_items + (m_start + i)) T(item);
		}

		void Insert(int i, T&& item) {
			if (m_start == 0) {
				if (m_count + m_start == m_capacity) {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, i);
					Move(newItems + i + 1, m_items + m_start + i, m_count - i);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
					m_count += 1;
				}
				else {
					Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
					m_count += 1;
				}
			}
			else if (m_count + m_start == m_capacity || i < m_count * 0.5) {
				Move(m_items + m_start + i - 1, m_items + m_start + i, i);
				m_start -= 1;
				m_count += 1;
			}
			else {
				Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
				m_count += 1;
			}
			new(m_items + (m_start + i)) T(std::forward<T>(item));
		}
		
		void Remove(int i) {
			FreeItem(m_items + (m_start + i));
			if (i < m_count * 0.5) {
				Move(m_items + m_start + 1, m_items + m_start, i);
				m_start += 1;
				m_count -= 1;
			}
			else {
				Move(m_items + m_start + i, m_items + m_start + i + 1, m_count - i - 1);
				m_count -= 1;
			}
		}

		void PopFirst() {
			FreeItem(m_items + m_start);
			m_start += 1;
			m_count -= 1;
		}

		void PopLast() {
			FreeItem(m_items + (m_start + m_count - 1));
			m_count -= 1;
		}

		int GetCount() const {
			return m_count;
		}

		T Get(int i) const {
			return m_items[i + m_start];
		}

		void Set(int i, const T& item) {
			T* p = m_items + (m_start + i);
			p->~T();
			new(p) T(item);
		}

		void Set(int i, T&& item) {
			T* p = m_items + (m_start + i);
			p->~T();
			new(p) T(std::forward<T>(item));
		}

		T* GetPointer(int i) {
			return m_items + (i + m_start);
		}

		const T* GetPointer(int i) const {
			return m_items + (i + m_start);
		}
	private:
		int CalculateCapacity(int count) {
			return count / 2 * 3 + 4;
		}

		void Move(T* dst, T* src, int count) {
			if (dst < src) {
				for (int i = 0; i < count; ++i) {
					new(dst + i) T(std::move(src[i]));
					(src + i)->~T();
				}
			}
			else if (dst > src) {
				for (int i = count - 1; i >= 0; --i) {
					new(dst + i) T(std::move(src[i]));
					(src + i)->~T();
				}
			}
		}

		void Copy(T* dst, T* src, int count) {
			for (int i = 0; i < count; ++i) {
				new(dst + i) T(src[i]);
			}
		}

		void FreeItem(T* p) {
			p->~T();
		}

		void FreeItems(T* p, int count) {
			for (int i = 0; i < count; ++i) {
				(p + i)->~T();
			}
		}
	private:
		int m_capacity;
		T* m_items;
		int m_start;
		int m_count;
	};

	/*
	template<class T>
	class SimpleArray {
	public:
		SimpleArray() {
			m_capacity = 0;
			m_items = nullptr;
			m_start = 0;
			m_count = 0;
		}

		SimpleArray(const SimpleArray& right) {
			m_capacity = right.m_count;
			m_start = 0;
			m_count = right.m_count;
			m_items = (T*)malloc(m_capacity * sizeof(T));
			Copy(m_items, right.m_items + right.m_start, m_count);
		}

		SimpleArray(SimpleArray&& right) {
			m_capacity = right.m_capacity;
			m_start = right.m_start;
			m_count = right.m_count;
			m_items = right.m_items;
			right.m_capacity = 0;
			right.m_start = 0;
			right.m_count = 0;
			right.m_items = nullptr;
		}

		virtual ~SimpleArray() {
			free(m_items);
		}

		void Append(const T& item) {
			if (m_start + m_count == m_capacity) {
				if (m_start > 0) {
					Move(m_items, m_items + m_start, m_count);
					m_start = 0;
				}
				else {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, m_count);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
				}
			}
			m_items[m_start + m_count] = item;
			m_count += 1;
		}

		void Insert(int i, const T& item) {
			if (m_start == 0) {
				if (m_count + m_start == m_capacity) {
					m_capacity = CalculateCapacity(m_capacity + 1);
					T* newItems = (T*)malloc(m_capacity * sizeof(T));
					Move(newItems, m_items + m_start, i);
					Move(newItems + i + 1, m_items + m_start + i, m_count - i);
					delete[] m_items;
					m_items = newItems;
					m_start = 0;
					m_count += 1;
				}
				else {
					Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
					m_count += 1;
				}
			}
			else if (m_count + m_start == m_capacity || i < m_count * 0.5) {
				Move(m_items + m_start + i - 1, m_items + m_start + i, i);
				m_start -= 1;
				m_count += 1;
			}
			else {
				Move(m_items + m_start + i + 1, m_items + m_start + i, m_count - i);
				m_count += 1;
			}
			m_items[m_start + i] = item;
		}

		void Remove(int i) {
			if (i < m_count * 0.5) {
				Move(m_items + m_start + 1, m_items + m_start, i);
				m_start += 1;
				m_count -= 1;
			}
			else {
				Move(m_items + m_start + i, m_items + m_start + i + 1, m_count - i - 1);
				m_count -= 1;
			}
		}

		void PopFirst() {
			m_start += 1;
			m_count -= 1;
		}

		void PopLast() {
			m_count -= 1;
		}

		int GetCount() const {
			return m_count;
		}

		T Get(int i) const {
			return m_items[i + m_start];
		}

		T* GetPointer(int i) {
			return m_items + (i + m_start);
		}

		const T* GetPointer(int i) const {
			return m_items + (i + m_start);
		}
	private:
		int CalculateCapacity(int count) {
			return count / 2 * 3 + 4;
		}

		void Move(T* dst, T* src, int count) {
			memcpy(dst, src, count * sizeof(T));
		}

		void Copy(T* dst, T* src, int count) {
			memcpy(dst, src, count * sizeof(T));
		}
	private:
		int m_capacity;
		T* m_items;
		int m_start;
		int m_count;
	};
	*/
}

#endif
