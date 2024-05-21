/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_STD_THREAD_
#define _WGP_STD_THREAD_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include "array.h"

namespace wgp {

    class Promise {
    public:
        Promise() : m_ready(false) {}
        virtual ~Promise() {}
    public:
        void SetReady() {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_ready = true;
            }
            m_condition.notify_all();
        }

        void Wait() {
            std::unique_lock<std::mutex> lock(m_mutex);
            while (!m_ready) {
                m_condition.wait(lock);
            }
        }
    private:
        bool m_ready;
        std::mutex m_mutex;
        std::condition_variable m_condition;
    };

    class Spinlock {
    public:
        void spinlock() {
            while (locked.test_and_set(std::memory_order_acquire)) {
            }
        }

        void lock() {
            while (locked.test_and_set(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
        }

        void unlock() {
            locked.clear(std::memory_order_release);
        }
    private:
        std::atomic_flag locked = ATOMIC_FLAG_INIT;
    };

    struct ThreadPoolTask {
        std::function<void()> Func;
        ThreadPoolTask* Next;
    };

    class ThreadPool {
    public:
        ThreadPool(int thread_count) : m_terminate(false) {
            m_thread_count = thread_count;
            m_threads = (std::thread*)malloc(m_thread_count * sizeof(std::thread));
            for (int i = 0; i < m_thread_count; ++i) {
                new(m_threads + i) std::thread(ThreadPool::ThreadFunc, this);
            }
            m_tail = nullptr;
            m_first_cache = nullptr;
        }

        virtual ~ThreadPool() {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_terminate = true;
            }
            m_condition.notify_all();
            for (int i = 0; i < m_thread_count; ++i) {
                m_threads[i].join();
            }
        }

        template<class F, class... Args>
        void AddTask(F&& f, Args&&... args) {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                ThreadPoolTask* task;
                if (m_first_cache) {
                    task = m_first_cache;
                    m_first_cache = task->Next;
                }
                else {
                    task = new ThreadPoolTask();
                }
                task->Func = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
                if (m_tail) {
                    task->Next = m_tail->Next;
                    m_tail->Next = task;
                    m_tail = task;
                }
                else {
                    task->Next = task;
                    m_tail = task;
                }
            }
            m_condition.notify_one();
        }
    private:
        static void ThreadFunc(ThreadPool* thread_pool) {
            while (true) {
                std::function<void()> func;
                {
                    std::unique_lock<std::mutex> lock(thread_pool->m_mutex);
                    while (!thread_pool->m_terminate && !thread_pool->m_tail) {
                        thread_pool->m_condition.wait(lock);
                    }
                    if (thread_pool->m_terminate && !thread_pool->m_tail) {
                        break;
                    }
                    ThreadPoolTask* task;
                    if (thread_pool->m_tail->Next == thread_pool->m_tail) {
                        task = thread_pool->m_tail;
                        thread_pool->m_tail = nullptr;
                    }
                    else {
                        task = thread_pool->m_tail->Next;
                        thread_pool->m_tail->Next = task->Next;
                    }
                    func = task->Func;
                    task->Next = thread_pool->m_first_cache;
                    thread_pool->m_first_cache = task;
                }
                func();
            }
        }
    private:
        ThreadPoolTask* m_tail;
        int m_thread_count;
        std::thread* m_threads;
        std::mutex m_mutex;
        std::condition_variable m_condition;
        bool m_terminate;
    private:
        ThreadPoolTask* m_first_cache;
    public:
        static ThreadPool* GetInstance() {
            return m_instance;
        }
    private:
        static ThreadPool* m_instance;
        class InstanceDestroyer {
        public:
            ~InstanceDestroyer();
        };
        static InstanceDestroyer m_instance_destroyer;
    };

    class ParallelComputor {
    public:
        ParallelComputor(int thread_count) : m_terminate(false) {
            m_state = 0;
            m_thread_count = thread_count;
            m_threads = (std::thread*)malloc(m_thread_count * sizeof(std::thread));
            m_current_indices = new int[m_thread_count];
            for (int i = 0; i < m_thread_count; ++i) {
                new(m_threads + i) std::thread(ParallelComputor::ThreadFunc, this, i);
                m_current_indices[i] = 0;
            }
        }

        virtual ~ParallelComputor() {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_terminate = true;
            }
            m_thread_condition.notify_all();
            for (int i = 0; i < m_thread_count; ++i) {
                m_threads[i].join();
            }
            delete[] m_current_indices;
        }

        void BeginAddTasks() {
            std::unique_lock<std::mutex> lock(m_mutex);
            while (m_state != 0) {
                m_add_condition.wait(lock);
            }
            m_state = -1;
        }

        template<class F, class... Args>
        void AddTask(F&& f, Args&&... args) {
            m_tasks.Append(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
        }

        void EndAddTasks() {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                if (m_tasks.GetCount() < m_thread_count) {
                    m_state = m_tasks.GetCount() + 1;
                }
                else {
                    m_state = m_thread_count + 1;
                }
                for (int i = 0; i < m_thread_count; ++i) {
                    m_current_indices[i] = i;
                }
            }
            m_thread_condition.notify_all();
        }

        void Wait() {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                while (m_state != 1) {
                    m_wait_condition.wait(lock);
                }
                m_state = 0;
                m_tasks.Clear();
            }
            m_add_condition.notify_all();
        }
    private:
        static void ThreadFunc(ParallelComputor* computor, int index) {
            while (true) {
                {
                    std::unique_lock<std::mutex> lock(computor->m_mutex);
                    while (!computor->m_terminate && (computor->m_state <= 1 || computor->m_current_indices[index] >= computor->m_tasks.GetCount())) {
                        computor->m_thread_condition.wait(lock);
                    }
                    if (computor->m_terminate && computor->m_current_indices[index] >= computor->m_tasks.GetCount()) {
                        break;
                    }
                }
                while (computor->m_current_indices[index] < computor->m_tasks.GetCount()) {
                    computor->m_tasks.Get(computor->m_current_indices[index])();
                    computor->m_current_indices[index] += computor->m_thread_count;
                }
                bool finish;
                {
                    std::unique_lock<std::mutex> lock(computor->m_mutex);
                    finish = --computor->m_state == 1;
                }
                if (finish) {
                    computor->m_wait_condition.notify_all();
                }
            }
        }
    private:
        int m_thread_count;
        std::thread* m_threads;
        int* m_current_indices;
        std::mutex m_mutex;
        std::condition_variable m_thread_condition;
        std::condition_variable m_add_condition;
        std::condition_variable m_wait_condition;
        int m_state;
        bool m_terminate;
        Array<std::function<void()>> m_tasks;
    public:
        static ParallelComputor* GetInstance() {
            return m_instance;
        }
    private:
        static ParallelComputor* m_instance;
        class InstanceDestroyer {
        public:
            ~InstanceDestroyer();
        };
        static InstanceDestroyer m_instance_destroyer;
    };

}

#endif