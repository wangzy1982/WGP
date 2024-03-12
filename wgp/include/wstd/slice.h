/*
	Original Author: Zuoyuan Wang
	Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_SLICE_
#define _WGP_STD_SLICE_

#include <atomic>
#include <memory>
#include "wbase.h"

struct SliceHeader {
	int Capacity;
	std::atomic<int> RefCount;
};

template<class T>
class Slice {
public:
	Slice();
	Slice(int capacity);
	Slice(const Slice& slice);
	Slice(const Slice& slice, int start, int count);
	virtual ~Slice();
	bool IsUnique() const;
	Slice operator=(const Slice& slice);
	int GetCount() const;
	void Resize(int start, int count);
	T GetItem(int index) const;
	T* GetItemPointer(int index) const;
	void SetItem(int index, const T& value);
	void Append(const T& value);
	void Append(const Slice& slice);
	void Insert(int index, const T& value);
	void Remove(int index);
	T Pop();
	void Clear();
private:
	void Decref();
	void Grow(int count);
private:
	SliceHeader* m_header;
	T* m_data;
	int m_start;
	int m_count;
};

template<class T> Slice<T>::Slice() {
	m_header = new SliceHeader();
	m_header->Capacity = 0;
	m_header->RefCount = 1;
	m_data = nullptr;
	m_start = 0;
	m_count = 0;
}

template<class T> Slice<T>::Slice(int capacity) {
	m_header = new SliceHeader();
	m_header->Capacity = capacity;
	m_header->RefCount = 1;
	m_data = new T[capacity];
	m_start = 0;
	m_count = 0;
}

template<class T> Slice<T>::Slice(const Slice& slice) {
	++slice.m_header->RefCount;
	m_header = slice.m_header;
	m_data = slice.m_data;
	m_start = slice.m_start;
	m_count = slice.m_count;
}

template<class T> Slice<T>::Slice(const Slice& slice, int start, int count) {
	++slice.m_header->RefCount;
	m_header = slice.m_header;
	m_data = slice.m_data;
	m_start = start;
	m_count = count;
}

template<class T> Slice<T>::~Slice() {
	Decref();
}

template<class T> bool Slice<T>::IsUnique() const {
	return m_header->RefCount == 1;
}

template<class T> Slice<T> Slice<T>::operator=(const Slice& slice) {
	++slice.m_header->RefCount;
	Decref();
	m_header = slice.m_header;
	m_data = slice.m_data;
	m_start = slice.m_start;
	m_count = slice.m_count;
	return *this;
}

template<class T> int Slice<T>::GetCount() const {
	return m_count;
}

template<class T> void Slice<T>::Resize(int start, int count) {
	m_start += start;
	m_count = count;
}

template<class T> T Slice<T>::GetItem(int index) const {
	return m_data[m_start + index];
}

template<class T> T* Slice<T>::GetItemPointer(int index) const {
	return &m_data[m_start + index];
}

template<class T> void Slice<T>::SetItem(int index, const T& value) {
	m_data[m_start + index] = value;
}

template<class T> void Slice<T>::Append(const T& value) {
	Grow(1);
	m_data[m_start + m_count] = value;
	++m_count;
}

template<class T> void Slice<T>::Insert(int index, const T& value) {
	Grow(1);
	for (int i = 0; i < m_count - index; ++i) {
		m_data[m_start + index + 1 + i] = m_data[m_start + index + i];
	}
	m_data[m_start + index] = value;
	++m_count;
}

template<class T> void Slice<T>::Remove(int index) {
	Grow(0);
	for (int i = 0; i < m_count - index - 1; ++i) {
		m_data[m_start + index + i] = m_data[m_start + index + 1 + i];
	}
	--m_count;
}

template<class T> void Slice<T>::Append(const Slice& slice) {
	if (slice.m_count > 0) {
		Grow(slice.m_count);
		for (int i = 0; i < slice.m_count; ++i) {
			m_data[m_start + m_count + i] = slice.m_data[slice.m_start + i];
		}
		m_count += slice.m_count;
	}
}

template<class T> T Slice<T>::Pop() {
	--m_count;
	return m_data[m_start + m_count];
}

template<class T> void Slice<T>::Clear() {
	m_count = 0;
}

template<class T> void Slice<T>::Decref() {
	if (--m_header->RefCount == 0) {
		delete[] m_data;
		delete m_header;
	}
}

template<class T> void Slice<T>::Grow(int count) {
	if (m_header->RefCount == 1) {
		int total_count = m_start + m_count + count;
		if (total_count > m_header->Capacity) {
			if (total_count < 16) {
				m_header->Capacity = total_count + 4;
			}
			else {
				m_header->Capacity = total_count + total_count / 4;
			}
			T* new_data = new T[m_header->Capacity];
			for (int i = 0; i < m_count; ++i) {
				new_data[i] = m_data[m_start + i];
			}
			delete[] m_data;
			m_data = new_data;
			m_start = 0;
		}
	}
	else {
		int new_capacity = m_header->Capacity;
		int total_count = m_start + m_count + count;
		if (total_count > new_capacity) {
			if (total_count < 16) {
				new_capacity = total_count + 4;
			}
			else {
				new_capacity = total_count + total_count / 4;
			}
		}
		T* new_data = new T[new_capacity];
		for (int i = 0; i < m_count; ++i) {
			new_data[i] = m_data[m_start + i];
		}
		Decref();
		m_header = new SliceHeader();
		m_header->Capacity = new_capacity;
		m_header->RefCount = 1;
		m_data = new_data;
		m_start = 0;
	}
}

template class WGP_API Slice<int>;
template class WGP_API Slice<double>;
template class WGP_API Slice<float>;

#endif
