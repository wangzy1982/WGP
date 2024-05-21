/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#include "wstd/thread.h"

namespace wgp {

    ThreadPool* ThreadPool::m_instance = new ThreadPool(8);

    ThreadPool::InstanceDestroyer::~InstanceDestroyer() {
        delete ThreadPool::GetInstance();
    };
    
    ThreadPool::InstanceDestroyer ThreadPool::m_instance_destroyer;

    ParallelComputor* ParallelComputor::m_instance = new ParallelComputor(8);

    ParallelComputor::InstanceDestroyer::~InstanceDestroyer() {
        delete ParallelComputor::GetInstance();
    };

    ParallelComputor::InstanceDestroyer ParallelComputor::m_instance_destroyer;
}