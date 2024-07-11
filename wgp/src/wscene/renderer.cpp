/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/

#include <assert.h>
#include "wscene/renderer.h"

namespace wgp {

    RenderingObject::RenderingObject(int classification, RenderingMaterial* material) : 
        m_classification(classification),
        m_material(material) {
    }

    RenderingMaterial* RenderingObject::GetMaterial() const {
        return m_material;
    }

    int RenderingObject::GetClassification() const {
        return m_classification;
    }

    RenderingObject::~RenderingObject() {
        delete m_material;
    }

    RenderingObjectFragment::RenderingObjectFragment(RenderingObject* rendering_object) :
        m_rendering_object(rendering_object) {
        if (m_rendering_object) {
            m_rendering_object->IncRef();
        }
    }

    RenderingObjectFragment::~RenderingObjectFragment() {
        if (m_rendering_object) {
            m_rendering_object->DecRef();
        }
    }

    RenderingObject* RenderingObjectFragment::GetRenderingObject() const {
        return m_rendering_object;
    }

    int RenderingObjectFragment::GetClassification() const {
        return m_rendering_object->GetClassification();
    }

    NullRenderingObjectFragment::NullRenderingObjectFragment(int classification) : 
        RenderingObjectFragment(nullptr),
        m_classification(classification) {
    }

    int NullRenderingObjectFragment::GetClassification() const {
        return m_classification;
    }

    void NullRenderingObjectFragment::SetState(int state) {
    }

    TYPE_IMP_0(RenderingMaterial);

    RenderingTree::RenderingTree(Model* model, bool is_order_affected_rendering, int complexity) :
        m_model(model),
        m_is_order_affected_rendering(is_order_affected_rendering),
        m_feature_info_root(nullptr),
        m_rendering_group_root(nullptr),
        m_dirty_level(2),
        m_complexity(complexity) {
        m_model->IncRef();
    }

    RenderingTree::~RenderingTree() {
        if (m_rendering_group_root) {
            FreeRenderingGroups(m_rendering_group_root);
        }
        if (m_feature_info_root) {
            FreeFeatureInfos(m_feature_info_root);
        }
        m_model->DecRef();
    }

    void RenderingTree::Notify(const Array<CommandLog*>& logs) {
        Array<Feature*> affected_features;
        for (int i = 0; i < logs.GetCount(); ++i) {
            CommandLog* log = logs.Get(i);
            int dirty_level = GetDirtyLevel(log);
            if (dirty_level > 0) {
                log->AppendAffectedFeature(affected_features);
                for (int j = 0; j < affected_features.GetCount(); ++j) {
                    SetFeatureInfoDirty(affected_features.Get(j));
                }
                affected_features.Clear();
                if (dirty_level == 2) {
                    m_dirty_level = 2;
                }
            }
        }
    }

    Model* RenderingTree::GetModel() const {
        return m_model;
    }

    void RenderingTree::RemoveRenderingObjects(int classification) {
        if (m_rendering_group_root) {
            RemoveRenderingObjects(m_rendering_group_root, classification);
        }
    }

    void RenderingTree::Render(Renderer* renderer, int classification) {
        Update();
        if (m_rendering_group_root) {
            Render(m_rendering_group_root, renderer, classification);
        }
    }

    void RenderingTree::SetState(Feature* feature, const Array<wgp::Feature*>& path, int state) {
        Update();
        FeatureInfo* feature_info = GetFeatureInfo(feature, path, false);
        if (feature_info) {
            for (int i = 0; i < feature_info->Fragments.GetCount(); ++i) {
                feature_info->Fragments.Get(i)->SetState(state);
            }
        }
    }

    void RenderingTree::ClearState(int clear_state) {
        Update();
        if (m_feature_info_root) {
            ClearState(m_feature_info_root, clear_state);
        }
    }

    void RenderingTree::Update() {
        if (m_dirty_level == 1) {
            UpdateDirty1(m_feature_info_root);
            if (m_rendering_group_root) {
                BuildRenderingObjects(m_rendering_group_root);
            }
            m_dirty_level = 0;
        }
        else if (m_dirty_level == 2) {
            Array<Feature*> path;
            Matrix4x4 matrix = Matrix4x4::Identity();
            int order = 1;
            UpdateDirty2(path, matrix, order, m_model);
            RemoveUnrenderedFeatureInfo(&m_feature_info_root);
            if (m_rendering_group_root) {
                BuildRenderingObjects(m_rendering_group_root);
            }
            m_dirty_level = 0;
        }
    }

    void RenderingTree::UpdateDirty1(FeatureInfo* feature_info) {
        if (feature_info->IsDirty) {
            UpdateDirtyFeatureInfo(feature_info);
        }
        if (feature_info->LeftChild) {
            UpdateDirty1(feature_info->LeftChild);
        }
        if (feature_info->RightChild) {
            UpdateDirty1(feature_info->RightChild);
        }
    }

    void RenderingTree::UpdateDirty2(Array<Feature*>& path, const Matrix4x4& matrix, int& order, Model* model) {
        Array<Feature*> features(model->GetFeatureCount());
        for (int i = 0; i < model->GetFeatureCount(); ++i) {
            Feature* feature = model->GetFeature(i);
            if (IsRenderingFeature(feature)) {
                features.Append(feature);
            }
        }
        SortFeatures(features);
        for (int i = 0; i < features.GetCount(); ++i) {
            Feature* feature = features.Get(i);
            FeatureInfo* feature_info = GetFeatureInfo(feature, path, true);
            if (order > feature_info->Order) {
                feature_info->Order = order;
                ++order;
                if (m_is_order_affected_rendering) {
                    feature_info->IsDirty = true;
                    for (int k = 0; k < feature_info->Fragments.GetCount(); ++k) {
                        delete feature_info->Fragments.Get(k);
                    }
                    feature_info->Fragments.Clear();
                }
            }
            if (feature_info->IsDirty) {
                UpdateDirtyFeatureInfo(feature_info);
                feature_info->IsDirty = false;
            }
            feature_info->IsRendered = true;
            Matrix4x4 temp_matrix;
            Model* model2 = GetReferenceModel(feature, temp_matrix);
            if (model2) {
                path.Append(feature);
                Matrix4x4 matrix2 = matrix * temp_matrix;
                UpdateDirty2(path, matrix2, order, model2);
                path.PopLast();
            }
        }
    }

    void RenderingTree::RemoveUnrenderedFeatureInfo(FeatureInfo** feature_info_address) {
        FeatureInfo* feature_info = *feature_info_address;
        if (!feature_info) {
            return;
        }
        if (feature_info->IsRendered) {
            feature_info->IsRendered = false;
            RemoveUnrenderedFeatureInfo(&feature_info->LeftChild);
            RemoveUnrenderedFeatureInfo(&feature_info->RightChild);
            return;
        }
        if (!feature_info->LeftChild) {
            *feature_info_address = feature_info->RightChild;
        }
        else if (!feature_info->RightChild) {
            *feature_info_address = feature_info->LeftChild;
        }
        else {
            FeatureInfo* feature_info2 = feature_info->LeftChild;
            if (!feature_info2->RightChild) {
                feature_info2->RightChild = feature_info->RightChild;
                *feature_info_address = feature_info2;
            }
            else {
                FeatureInfo** feature_info_address2 = &feature_info2->RightChild;
                while ((*feature_info_address2)->RightChild) {
                    feature_info_address2 = &(*feature_info_address2)->RightChild;
                }
                FeatureInfo* feature_info2 = *feature_info_address2;
                *feature_info_address2 = feature_info2->LeftChild;
                feature_info2->LeftChild = feature_info->LeftChild;
                feature_info2->RightChild = feature_info->RightChild;
                *feature_info_address = feature_info2;
            }
        }
        RemoveFeatureInfo(feature_info);
    }

    void RenderingTree::UpdateDirtyFeatureInfo(FeatureInfo* feature_info) {
        RemoveFeatureInfoFromGroup(feature_info);
        RenderingMaterial* material = nullptr;
        bool is_classification_enabled = false;
        Calculate(feature_info, material, is_classification_enabled, feature_info->Box, feature_info->Complexity);
        if (material) {
            feature_info->Group = GetRenderingGroup(material, is_classification_enabled, true);
            AddFeatureInfoToGroup(feature_info);
        }
    }

    void RenderingTree::RemoveFeatureInfo(FeatureInfo* feature_info) {
        RemoveFeatureInfoFromGroup(feature_info);
        for (int i = 0; i < feature_info->Fragments.GetCount(); ++i) {
            delete feature_info->Fragments.Get(i);
        }
        for (int i = 0; i < feature_info->Path.GetCount(); ++i) {
            feature_info->Path.Get(i)->DecRef();
        }
        feature_info->Feature->DecRef();
        delete feature_info;
    }

    void RenderingTree::AddFeatureInfoToGroup(FeatureInfo* feature_info) {
        if (!feature_info->Group) {
            return;
        }
        if (!feature_info->Group->Root) {
            RenderingNode* new_root = new RenderingNode;
            new_root->Height = 0;
            new_root->Box = feature_info->Box;
            new_root->Complexity = feature_info->Complexity;
            new_root->Children[0] = feature_info;
            new_root->Children[1] = nullptr;
            new_root->Children[2] = nullptr;
            new_root->IsDirty = true;
            feature_info->Group->Root = new_root;
            return;
        }
        RenderingNode* new_node = AddFeatureInfoToGroup(feature_info->Group->Root, feature_info);
        if (new_node) {
            RenderingNode* new_root = new RenderingNode;
            new_root->Height = new_node->Height + 1;
            new_root->Box = new_node->Box;
            new_root->Box.Merge(feature_info->Group->Root->Box);
            new_root->Complexity = new_node->Complexity + feature_info->Group->Root->Complexity;
            new_root->Children[0] = feature_info->Group->Root;
            new_root->Children[1] = new_node;
            new_root->Children[2] = nullptr;
            new_root->IsDirty = true;
            feature_info->Group->Root = new_root;
        }
    }

    RenderingNode* RenderingTree::AddFeatureInfoToGroup(RenderingNode* node, FeatureInfo* feature_info) {
        if (node->Height == 0) {
            if (node->Children[1] == nullptr) {
                node->Children[1] = feature_info;
                node->Box.Merge(feature_info->Box);
                node->Complexity += feature_info->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            if (node->Children[2] == nullptr) {
                node->Children[2] = feature_info;
                node->Box.Merge(feature_info->Box);
                node->Complexity += feature_info->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            FeatureInfo* feature_infos[4] = {
                feature_info,
                (FeatureInfo*)node->Children[0],
                (FeatureInfo*)node->Children[1],
                (FeatureInfo*)node->Children[2]
            };
            int mi = -1;
            int mj = -1;
            double md = -1;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    Interval3d box = feature_infos[i]->Box;
                    box.Merge(feature_infos[j]->Box);
                    double d = box.DiagonalLength();
                    if (d > md) {
                        md = d;
                        mi = i;
                        mj = j;
                    }
                }
            }
            if (mi != 0) {
                FeatureInfo* t = feature_infos[0];
                feature_infos[0] = feature_infos[mi];
                feature_infos[mi] = t;
            }
            if (mj != 1) {
                FeatureInfo* t = feature_infos[1];
                feature_infos[1] = feature_infos[mj];
                feature_infos[mj] = t;
            }
            Interval3d box = feature_infos[0]->Box;
            box.Merge(feature_infos[2]->Box);
            double d1 = box.DiagonalLength();
            box = feature_infos[1]->Box;
            box.Merge(feature_infos[3]->Box);
            d1 += box.DiagonalLength();
            box = feature_infos[0]->Box;
            box.Merge(feature_infos[3]->Box);
            double d2 = box.DiagonalLength();
            box = feature_infos[1]->Box;
            box.Merge(feature_infos[2]->Box);
            d2 += box.DiagonalLength();
            if (d1 > d2) {
                FeatureInfo* t = feature_infos[2];
                feature_infos[2] = feature_infos[3];
                feature_infos[3] = t;
            }
            node->Box = feature_infos[0]->Box;
            node->Box.Merge(feature_infos[2]->Box);
            node->Complexity = feature_infos[0]->Complexity + feature_infos[2]->Complexity;
            node->Children[0] = feature_infos[0];
            node->Children[1] = feature_infos[2];
            node->Children[2] = nullptr;
            node->IsDirty = true;
            for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                node->RenderingObjects.Get(k)->DecRef();
            }
            node->RenderingObjects.Clear();
            RenderingNode* new_node = new RenderingNode;
            new_node->Height = node->Height;
            new_node->Box = feature_infos[1]->Box;
            new_node->Box.Merge(feature_infos[3]->Box);
            new_node->Complexity = feature_infos[1]->Complexity + feature_infos[3]->Complexity;
            new_node->Children[0] = feature_infos[1];
            new_node->Children[1] = feature_infos[3];
            new_node->Children[2] = nullptr;
            new_node->IsDirty = true;
            return new_node;
        }
        else {
            int mi = -1;
            double md = 1E100;
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                Interval3d box = ((RenderingNode*)node->Children[i])->Box;
                box.Merge(feature_info->Box);
                double d = box.DiagonalLength();
                if (d < md) {
                    md = d;
                    mi = i;
                }
            }
            assert(mi >= 0 && mi < 3);
            RenderingNode* new_child_node = AddFeatureInfoToGroup((RenderingNode*)node->Children[mi], feature_info);
            if (!new_child_node) {
                node->Box.Merge(feature_info->Box);
                node->Complexity += feature_info->Complexity;
                return nullptr;
            }
            if (node->Children[1] == nullptr) {
                node->Children[1] = new_child_node;
                node->Box.Merge(feature_info->Box);
                node->Complexity += feature_info->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            if (node->Children[2] == nullptr) {
                node->Children[2] = new_child_node;
                node->Box.Merge(feature_info->Box);
                node->Complexity += feature_info->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            RenderingNode* nodes[4] = {
                new_child_node,
                (RenderingNode*)node->Children[0],
                (RenderingNode*)node->Children[1],
                (RenderingNode*)node->Children[2]
            };
            mi = -1;
            int mj = -1;
            md = -1;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    Interval3d box = nodes[i]->Box;
                    box.Merge(nodes[j]->Box);
                    double d = box.DiagonalLength();
                    if (d > md) {
                        md = d;
                        mi = i;
                        mj = j;
                    }
                }
            }
            if (mi != 0) {
                RenderingNode* t = nodes[0];
                nodes[0] = nodes[mi];
                nodes[mi] = t;
            }
            if (mj != 1) {
                RenderingNode* t = nodes[1];
                nodes[1] = nodes[mj];
                nodes[mj] = t;
            }
            Interval3d box = nodes[0]->Box;
            box.Merge(nodes[2]->Box);
            double d1 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[3]->Box);
            d1 += box.DiagonalLength();
            box = nodes[0]->Box;
            box.Merge(nodes[3]->Box);
            double d2 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[2]->Box);
            d2 += box.DiagonalLength();
            if (d1 > d2) {
                RenderingNode* t = nodes[2];
                nodes[2] = nodes[3];
                nodes[3] = t;
            }
            node->Box = nodes[0]->Box;
            node->Box.Merge(nodes[2]->Box);
            node->Complexity = nodes[0]->Complexity + nodes[2]->Complexity;
            node->Children[0] = nodes[0];
            node->Children[1] = nodes[2];
            node->Children[2] = nullptr;
            node->IsDirty = true;
            for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                node->RenderingObjects.Get(k)->DecRef();
            }
            node->RenderingObjects.Clear();
            RenderingNode* new_node = new RenderingNode;
            new_node->Height = node->Height;
            new_node->Box = nodes[1]->Box;
            new_node->Box.Merge(nodes[3]->Box);
            new_node->Complexity = nodes[1]->Complexity + nodes[3]->Complexity;
            new_node->Children[0] = nodes[1];
            new_node->Children[1] = nodes[3];
            new_node->Children[2] = nullptr;
            new_node->IsDirty = true;
            return new_node;
        }
    }

    RenderingNode* RenderingTree::ReinsertRenderingNodeToGroup(RenderingNode* node, RenderingNode* reinsert_node) {
        if (node->Height == reinsert_node->Height + 1) {
            if (node->Children[1] == nullptr) {
                node->Children[1] = reinsert_node;
                node->Box.Merge(reinsert_node->Box);
                node->Complexity += reinsert_node->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            if (node->Children[2] == nullptr) {
                node->Children[2] = reinsert_node;
                node->Box.Merge(reinsert_node->Box);
                node->Complexity += reinsert_node->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            RenderingNode* nodes[4] = {
                reinsert_node,
                (RenderingNode*)node->Children[0],
                (RenderingNode*)node->Children[1],
                (RenderingNode*)node->Children[2]
            };
            int mi = -1;
            int mj = -1;
            double md = -1;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    Interval3d box = nodes[i]->Box;
                    box.Merge(nodes[j]->Box);
                    double d = box.DiagonalLength();
                    if (d > md) {
                        md = d;
                        mi = i;
                        mj = j;
                    }
                }
            }
            if (mi != 0) {
                RenderingNode* t = nodes[0];
                nodes[0] = nodes[mi];
                nodes[mi] = t;
            }
            if (mj != 1) {
                RenderingNode* t = nodes[1];
                nodes[1] = nodes[mj];
                nodes[mj] = t;
            }
            Interval3d box = nodes[0]->Box;
            box.Merge(nodes[2]->Box);
            double d1 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[3]->Box);
            d1 += box.DiagonalLength();
            box = nodes[0]->Box;
            box.Merge(nodes[3]->Box);
            double d2 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[2]->Box);
            d2 += box.DiagonalLength();
            if (d1 > d2) {
                RenderingNode* t = nodes[2];
                nodes[2] = nodes[3];
                nodes[3] = t;
            }
            node->Box = nodes[0]->Box;
            node->Box.Merge(nodes[2]->Box);
            node->Complexity = nodes[0]->Complexity + nodes[2]->Complexity;
            node->Children[0] = nodes[0];
            node->Children[1] = nodes[2];
            node->Children[2] = nullptr;
            node->IsDirty = true;
            for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                node->RenderingObjects.Get(k)->DecRef();
            }
            node->RenderingObjects.Clear();
            RenderingNode* new_node = new RenderingNode;
            new_node->Height = node->Height;
            new_node->Box = nodes[1]->Box;
            new_node->Box.Merge(nodes[3]->Box);
            new_node->Complexity = nodes[1]->Complexity + nodes[3]->Complexity;
            new_node->Children[0] = nodes[1];
            new_node->Children[1] = nodes[3];
            new_node->Children[2] = nullptr;
            new_node->IsDirty = true;
            return new_node;
        }
        else {
            int mi = -1;
            double md = 1E100;
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                Interval3d box = ((RenderingNode*)node->Children[i])->Box;
                box.Merge(reinsert_node->Box);
                double d = box.DiagonalLength();
                if (d < md) {
                    md = d;
                    mi = i;
                }
            }
            assert(mi >= 0 && mi < 3);
            RenderingNode* new_child_node = ReinsertRenderingNodeToGroup((RenderingNode*)node->Children[mi], reinsert_node);
            if (!new_child_node) {
                node->Box.Merge(reinsert_node->Box);
                node->Complexity += reinsert_node->Complexity;
                return nullptr;
            }
            if (node->Children[1] == nullptr) {
                node->Children[1] = new_child_node;
                node->Box.Merge(reinsert_node->Box);
                node->Complexity += reinsert_node->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            if (node->Children[2] == nullptr) {
                node->Children[2] = new_child_node;
                node->Box.Merge(reinsert_node->Box);
                node->Complexity += reinsert_node->Complexity;
                node->IsDirty = true;
                for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                    node->RenderingObjects.Get(k)->DecRef();
                }
                node->RenderingObjects.Clear();
                return nullptr;
            }
            RenderingNode* nodes[4] = {
                new_child_node,
                (RenderingNode*)node->Children[0],
                (RenderingNode*)node->Children[1],
                (RenderingNode*)node->Children[2]
            };
            mi = -1;
            int mj = -1;
            md = -1;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 4; ++j) {
                    Interval3d box = nodes[i]->Box;
                    box.Merge(nodes[j]->Box);
                    double d = box.DiagonalLength();
                    if (d > md) {
                        md = d;
                        mi = i;
                        mj = j;
                    }
                }
            }
            if (mi != 0) {
                RenderingNode* t = nodes[0];
                nodes[0] = nodes[mi];
                nodes[mi] = t;
            }
            if (mj != 1) {
                RenderingNode* t = nodes[1];
                nodes[1] = nodes[mj];
                nodes[mj] = t;
            }
            Interval3d box = nodes[0]->Box;
            box.Merge(nodes[2]->Box);
            double d1 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[3]->Box);
            d1 += box.DiagonalLength();
            box = nodes[0]->Box;
            box.Merge(nodes[3]->Box);
            double d2 = box.DiagonalLength();
            box = nodes[1]->Box;
            box.Merge(nodes[2]->Box);
            d2 += box.DiagonalLength();
            if (d1 > d2) {
                RenderingNode* t = nodes[2];
                nodes[2] = nodes[3];
                nodes[3] = t;
            }
            node->Box = nodes[0]->Box;
            node->Box.Merge(nodes[2]->Box);
            node->Complexity = nodes[0]->Complexity + nodes[2]->Complexity;
            node->Children[0] = nodes[0];
            node->Children[1] = nodes[2];
            node->Children[2] = nullptr;
            node->IsDirty = true;
            for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                node->RenderingObjects.Get(k)->DecRef();
            }
            node->RenderingObjects.Clear();
            RenderingNode* new_node = new RenderingNode;
            new_node->Height = node->Height;
            new_node->Box = nodes[1]->Box;
            new_node->Box.Merge(nodes[3]->Box);
            new_node->Complexity = nodes[1]->Complexity + nodes[3]->Complexity;
            new_node->Children[0] = nodes[1];
            new_node->Children[1] = nodes[3];
            new_node->Children[2] = nullptr;
            new_node->IsDirty = true;
            return new_node;
        }
    }

    bool RenderingTree::RemoveFeatureInfoFromGroup(FeatureInfo* feature_info) {
        if (!feature_info->Group) {
            return false;
        }
        if (!feature_info->Group->Root) {
            return false;
        }
        RenderingNode* node = feature_info->Group->Root;
        if (!feature_info->Box.IsInner(node->Box)) {
            return false;
        }
        if (node->Height == 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    return false;
                }
                if (node->Children[i] == feature_info) {
                    int j = i + 1;
                    while (j < 3) {
                        if (!node->Children[j]) {
                            break;
                        }
                        ++j;
                    }
                    --j;
                    if (i != j) {
                        node->Children[i] = node->Children[j];
                    }
                    node->Children[j] = nullptr;
                    if (j == 0) {
                        for (int i = 0; i < node->RenderingObjects.GetCount(); ++i) {
                            delete node->RenderingObjects.Get(i);
                        }
                        delete node;
                        feature_info->Group->Root = nullptr;
                    }
                    else {
                        node->Box = ((FeatureInfo*)node->Children[0])->Box;
                        node->Complexity = ((FeatureInfo*)node->Children[0])->Complexity;
                        if (node->Children[1]) {
                            node->Box.Merge(((FeatureInfo*)node->Children[1])->Box);
                            node->Complexity += ((FeatureInfo*)node->Children[1])->Complexity;
                        }
                        node->IsDirty = true;
                        for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                            node->RenderingObjects.Get(k)->DecRef();
                        }
                        node->RenderingObjects.Clear();
                    }
                    return true;
                }
            }
        }
        else {
            Array<RenderingNode*> reinsert_nodes;
            FeatureInfo* reinsert_feature_info = nullptr;
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    return false;
                }
                RenderingNode* child_node = (RenderingNode*)node->Children[i];
                if (RemoveFeatureInfoFromGroup(child_node, feature_info, reinsert_nodes, reinsert_feature_info)) {
                    if (child_node->Children[0]) {
                        node->Box = ((RenderingNode*)node->Children[0])->Box;
                        node->Complexity = ((RenderingNode*)node->Children[0])->Complexity;
                        for (int j = 1; j < 3; ++j) {
                            if (!node->Children[j]) {
                                break;
                            }
                            node->Box.Merge(((RenderingNode*)node->Children[j])->Box);
                            node->Complexity += ((RenderingNode*)node->Children[j])->Complexity;
                        }
                        node->IsDirty = true;
                        for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                            node->RenderingObjects.Get(k)->DecRef();
                        }
                        node->RenderingObjects.Clear();
                    }
                    else {
                        int j = i + 1;
                        while (j < 3) {
                            if (!node->Children[j]) {
                                break;
                            }
                            ++j;
                        }
                        --j;
                        for (int k = 0; k < child_node->RenderingObjects.GetCount(); ++k) {
                            child_node->RenderingObjects.Get(k)->DecRef();
                        }
                        delete child_node;
                        if (i != j) {
                            node->Children[i] = node->Children[j];
                        }
                        node->Children[j] = nullptr;
                        assert(j != 0);
                        node->Box = ((FeatureInfo*)node->Children[0])->Box;
                        node->Complexity = ((FeatureInfo*)node->Children[0])->Complexity;
                        if (node->Children[1]) {
                            node->Box.Merge(((FeatureInfo*)node->Children[1])->Box);
                            node->Complexity += ((FeatureInfo*)node->Children[1])->Complexity;
                        }
                        node->IsDirty = true;
                        for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                            node->RenderingObjects.Get(k)->DecRef();
                        }
                        node->RenderingObjects.Clear();
                    }
                    for (int j = reinsert_nodes.GetCount() - 1; j >= 0; --j) {
                        ReinsertRenderingNodeToGroup(node, reinsert_nodes.Get(j));
                    }
                    if (reinsert_feature_info) {
                        AddFeatureInfoToGroup(node, reinsert_feature_info);
                    }
                    while (node->Height > 0 && node->Children[1] == nullptr) {
                        RenderingNode* temp = (RenderingNode*)node->Children[0];
                        assert(node->RenderingObjects.GetCount() == 0);
                        delete node;
                        node = temp;
                    }
                    feature_info->Group->Root = node;
                    return true;
                }
            }
        }
        return false;
    }

    bool RenderingTree::RemoveFeatureInfoFromGroup(RenderingNode* node, FeatureInfo* feature_info,
        Array<RenderingNode*>& reinsert_nodes, FeatureInfo*& reinsert_feature_info) {
        if (!feature_info->Box.IsInner(node->Box)) {
            return false;
        }
        if (node->Height == 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    return false;
                }
                if (node->Children[i] == feature_info) {
                    int j = i + 1;
                    while (j < 3) {
                        if (!node->Children[j]) {
                            break;
                        }
                        ++j;
                    }
                    --j;
                    if (i != j) {
                        node->Children[i] = node->Children[j];
                    }
                    node->Children[j] = nullptr;
                    assert(j > 0);
                    if (j == 1) {
                        reinsert_feature_info = (FeatureInfo*)node->Children[0];
                        node->Children[0] = nullptr;
                    }
                    else {
                        node->Box = ((FeatureInfo*)node->Children[0])->Box;
                        node->Complexity = ((FeatureInfo*)node->Children[0])->Complexity;
                        node->Box.Merge(((FeatureInfo*)node->Children[1])->Box);
                        node->Complexity += ((FeatureInfo*)node->Children[1])->Complexity;
                        node->IsDirty = true;
                        for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                            node->RenderingObjects.Get(k)->DecRef();
                        }
                        node->RenderingObjects.Clear();
                    }
                    return true;
                }
            }
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    return false;
                }
                RenderingNode* child_node = (RenderingNode*)node->Children[i];
                if (RemoveFeatureInfoFromGroup(child_node, feature_info, reinsert_nodes, reinsert_feature_info)) {
                    if (child_node->Children[0]) {
                        node->Box = ((RenderingNode*)node->Children[0])->Box;
                        node->Complexity = ((RenderingNode*)node->Children[0])->Complexity;
                        for (int j = 1; j < 3; ++j) {
                            if (!node->Children[j]) {
                                break;
                            }
                            node->Box.Merge(((RenderingNode*)node->Children[j])->Box);
                            node->Complexity += ((RenderingNode*)node->Children[j])->Complexity;
                        }
                        node->IsDirty = true;
                        for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                            node->RenderingObjects.Get(k)->DecRef();
                        }
                        node->RenderingObjects.Clear();
                    }
                    else {
                        int j = i + 1;
                        while (j < 3) {
                            if (!node->Children[j]) {
                                break;
                            }
                            ++j;
                        }
                        --j;
                        for (int k = 0; k < child_node->RenderingObjects.GetCount(); ++k) {
                            child_node->RenderingObjects.Get(k)->DecRef();
                        }
                        delete child_node;
                        if (i != j) {
                            node->Children[i] = node->Children[j];
                        }
                        node->Children[j] = nullptr;
                        assert(j != 0);
                        if (j == 1) {
                            reinsert_nodes.Append((RenderingNode*)node->Children[0]);
                            node->Children[0] = nullptr;
                        }
                        else {
                            node->Box = ((FeatureInfo*)node->Children[0])->Box;
                            node->Complexity = ((FeatureInfo*)node->Children[0])->Complexity;
                            node->Box.Merge(((FeatureInfo*)node->Children[1])->Box);
                            node->Complexity += ((FeatureInfo*)node->Children[1])->Complexity;
                            node->IsDirty = true;
                            for (int k = 0; k < node->RenderingObjects.GetCount(); ++k) {
                                node->RenderingObjects.Get(k)->DecRef();
                            }
                            node->RenderingObjects.Clear();
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }

    void RenderingTree::BuildRenderingObjects(RenderingGroup* group) {
        if (group->LeftChild) {
            BuildRenderingObjects(group->LeftChild);
        }
        if (group->RightChild) {
            BuildRenderingObjects(group->RightChild);
        }
        if (group->Classifications.GetCount() == 0) {
            return;
        }
        if (!group->Root) {
            return;
        }
        BuildRenderingObjects(group, group->Root);
    }

    void RenderingTree::BuildRenderingObjects(RenderingGroup* group, RenderingNode* node) {
        if (!node->IsDirty) {
            return;
        }
        if (node->Complexity <= m_complexity || node->Height == 0) {
            BuildRenderingObjects(group, node, node->RenderingObjects);
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                BuildRenderingObjects(group, (RenderingNode*)node->Children[i]);
            }
        }
        node->IsDirty = false;
    }

    void RenderingTree::BuildRenderingObjects(RenderingGroup* group, RenderingNode* node, Array<RenderingObject*>& rendering_objects) {
        if (node->Height == 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                for (int j = 0; j < group->Classifications.GetCount(); ++j) {
                    int classification;
                    if (group->IsClassificationEnabled) {
                        classification = group->Classifications.Get(j);
                    }
                    else {
                        classification = 0;
                    }
                    BuildRenderingObjects(classification, (FeatureInfo*)node->Children[i], rendering_objects);
                }
            }
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                RenderingNode* child_node = (RenderingNode*)node->Children[i];
                child_node->IsDirty = true;
                for (int k = 0; k < child_node->RenderingObjects.GetCount(); ++k) {
                    child_node->RenderingObjects.Get(k)->DecRef();
                }
                child_node->RenderingObjects.Clear();
                BuildRenderingObjects(group, child_node, rendering_objects);
            }
        }
    }

    void RenderingTree::BuildRenderingObjects(int classification, FeatureInfo* feature_info, Array<RenderingObject*>& rendering_objects) {        
        Array<RenderingObjectFragment*> fragments = std::move(feature_info->Fragments);
        feature_info->Fragments = Array<RenderingObjectFragment*>();
        bool is_new_fragment = true;
        for (int j = 0; j < fragments.GetCount(); ++j) {
            RenderingObjectFragment* fragment = fragments.Get(j);
            if (fragment->GetClassification() != classification) {
                continue;
            }
            if (!fragment->GetRenderingObject()) {
                feature_info->Fragments.Append(fragment);
                continue;
            }
            RenderingMaterial* material = fragment->GetRenderingObject()->GetMaterial();
            RenderingObject* rendering_object = nullptr;
            for (int k = 0; k < rendering_objects.GetCount(); ++k) {
                RenderingObject* rendering_object2 = rendering_objects.Get(k);
                if (rendering_object2->GetClassification() == classification) {
                    RenderingMaterial* material2 = rendering_object2->GetMaterial();
                    if (material == material2) {
                        rendering_object = rendering_object2;
                        break;
                    }
                    if (material && material2 && material->Compare(material2) == 0) {
                        rendering_object = rendering_object2;
                        break;
                    }
                }
            }
            if (rendering_object) {
                feature_info->Fragments.Append(rendering_object->Merge(fragment));
                delete fragment;
            }
            else {
                if (fragment->GetRenderingObject()) {
                    rendering_object = fragment->GetRenderingObject()->NewRenderingObject(fragment);
                    rendering_object->IncRef();
                    rendering_objects.Append(rendering_object);
                    feature_info->Fragments.Append(rendering_object->NewFragment());
                }
                else {
                    feature_info->Fragments.Append(new NullRenderingObjectFragment(classification));
                }
                delete fragment;
            }
            is_new_fragment = false;
        }
        if (is_new_fragment) {
            Array<RenderingObject*> current_rendering_objects;
            BuildRenderingObject(feature_info, classification, current_rendering_objects);
            if (current_rendering_objects.GetCount() == 0) {
                feature_info->Fragments.Append(new NullRenderingObjectFragment(classification));
            }
            else {
                for (int j = 0; j < current_rendering_objects.GetCount(); ++j) {
                    RenderingObject* current_rendering_object = current_rendering_objects.Get(j);
                    RenderingMaterial* material = current_rendering_object->GetMaterial();
                    RenderingObject* rendering_object = nullptr;
                    for (int k = 0; k < rendering_objects.GetCount(); ++k) {
                        RenderingObject* rendering_object2 = rendering_objects.Get(k);
                        if (rendering_object2->GetClassification() == classification) {
                            RenderingMaterial* material2 = rendering_object2->GetMaterial();
                            if (material == material2) {
                                rendering_object = rendering_object2;
                                break;
                            }
                            if (material && material2 && material->Compare(material2) == 0) {
                                rendering_object = rendering_object2;
                                break;
                            }
                        }
                    }
                    if (rendering_object) {
                        feature_info->Fragments.Append(rendering_object->Merge(current_rendering_object));
                        delete current_rendering_object;
                    }
                    else {
                        rendering_object = current_rendering_object;
                        rendering_object->IncRef();
                        rendering_objects.Append(rendering_object);
                        feature_info->Fragments.Append(rendering_object->NewFragment());
                    }
                }
            }
        }
    }

    void RenderingTree::AppendRenderingObjects(RenderingGroup* group, int classification) {
        if (group->Classifications.GetCount() == 0) {
            group->Classifications.Append(classification);
            BuildRenderingObjects(group);
            return;
        }
        if (!group->IsClassificationEnabled) {
            return;
        }
        if (!group->Root) {
            return;
        }
        for (int i = 0; i < group->Classifications.GetCount(); ++i) {
            if (group->Classifications.Get(i) == classification) {
                return;
            }
        }
        group->Classifications.Append(classification);
        AppendRenderingObjects(group, classification, group->Root);
    }

    void RenderingTree::AppendRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node) {
        if (node->Complexity <= m_complexity || node->Height == 0) {
            AppendRenderingObjects(group, classification, node, node->RenderingObjects);
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                AppendRenderingObjects(group, classification, (RenderingNode*)node->Children[i]);
            }
        }
    }

    void RenderingTree::AppendRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node, Array<RenderingObject*>& rendering_objects) {
        if (node->Height == 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                AppendRenderingObjects(classification, (FeatureInfo*)node->Children[i], rendering_objects);
            }
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                AppendRenderingObjects(group, classification, (RenderingNode*)node->Children[i], rendering_objects);
            }
        }
    }

    void RenderingTree::AppendRenderingObjects(int classification, FeatureInfo* feature_info, Array<RenderingObject*>& rendering_objects) {
        Array<RenderingObject*> current_rendering_objects;
        BuildRenderingObject(feature_info, classification, current_rendering_objects);
        if (current_rendering_objects.GetCount() == 0) {
            feature_info->Fragments.Append(new NullRenderingObjectFragment(classification));
        }
        else {
            for (int j = 0; j < current_rendering_objects.GetCount(); ++j) {
                RenderingObject* current_rendering_object = current_rendering_objects.Get(j);
                RenderingMaterial* material = current_rendering_object->GetMaterial();
                RenderingObject* rendering_object = nullptr;
                for (int k = 0; k < rendering_objects.GetCount(); ++k) {
                    RenderingObject* rendering_object2 = rendering_objects.Get(k);
                    if (rendering_object2->GetClassification() == classification) {
                        RenderingMaterial* material2 = rendering_object2->GetMaterial();
                        if (material == material2) {
                            rendering_object = rendering_object2;
                            break;
                        }
                        if (material && material2 && material->Compare(material2) == 0) {
                            rendering_object = rendering_object2;
                            break;
                        }
                    }
                }
                if (rendering_object) {
                    feature_info->Fragments.Append(rendering_object->Merge(current_rendering_object));
                    delete current_rendering_object;
                }
                else {
                    rendering_object = current_rendering_object;
                    rendering_object->IncRef();
                    rendering_objects.Append(rendering_object);
                    feature_info->Fragments.Append(rendering_object->NewFragment());
                }
            }
        }
    }

    void RenderingTree::RemoveRenderingObjects(RenderingGroup* group, int classification) {
        if (group->LeftChild) {
            RemoveRenderingObjects(group->LeftChild, classification);
        }
        if (group->RightChild) {
            RemoveRenderingObjects(group->RightChild, classification);
        }
        bool b = true;
        for (int i = 0; i < group->Classifications.GetCount(); ++i) {
            if (group->Classifications.Get(i) == classification) {
                group->Classifications.Remove(i);
                b = false;
                break;
            }
        }
        if (b) {
            return;
        }
        if (!group->IsClassificationEnabled && group->Classifications.GetCount() > 0) {
            return;
        }
        if (!group->Root) {
            return;
        }
        RemoveRenderingObjects(group, group->IsClassificationEnabled ? classification : 0, group->Root);
    }

    void RenderingTree::RemoveRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node) {
        int ci = 0;
        while (ci < node->RenderingObjects.GetCount()) {
            RenderingObject* rendering_object = node->RenderingObjects.Get(ci);
            if (rendering_object->GetClassification() == classification) {
                rendering_object->DecRef();
                if (ci != node->RenderingObjects.GetCount() - 1) {
                    node->RenderingObjects.Set(ci, node->RenderingObjects.Get(node->RenderingObjects.GetCount() - 1));
                }
                node->RenderingObjects.PopLast();
            }
            ++ci;
        }
        node->IsDirty = group->Classifications.GetCount() == 0;
        if (node->Height == 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                RemoveRenderingObjects(classification, (FeatureInfo*)node->Children[i]);
            }
        }
        else {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                RemoveRenderingObjects(group, classification, (RenderingNode*)node->Children[i]);
            }
        }
    }

    void RenderingTree::RemoveRenderingObjects(int classification, FeatureInfo* feature_info) {
        int ci = 0;
        while (ci < feature_info->Fragments.GetCount()) {
            RenderingObjectFragment* fragment = feature_info->Fragments.Get(ci);
            if (fragment->GetClassification() == classification) {
                delete fragment;
                if (ci != feature_info->Fragments.GetCount() - 1) {
                    feature_info->Fragments.Set(ci, feature_info->Fragments.Get(feature_info->Fragments.GetCount() - 1));
                }
                feature_info->Fragments.PopLast();
            }
            ++ci;
        }
    }

    void RenderingTree::Render(RenderingGroup* group, Renderer* renderer, int classification) {
        if (group->LeftChild) {
            Render(group->LeftChild, renderer, classification);
        }
        AppendRenderingObjects(group, classification);
        if (group->Root) {
            Render(group, group->Root, renderer, group->IsClassificationEnabled ? classification : 0);
        }
        if (group->RightChild) {
            Render(group->RightChild, renderer, classification);
        }
    }

    void RenderingTree::Render(RenderingGroup* group, RenderingNode* node, Renderer* renderer, int classification) {
        if (node->IsDirty) {
            return;
        }
        for (int i = 0; i < node->RenderingObjects.GetCount(); ++i) {
            RenderingObject* rendering_object = node->RenderingObjects.Get(i);
            if (rendering_object->GetClassification() == classification) {
                RenderingMaterial* material = rendering_object->GetMaterial();
                if (!material) {
                    material = group->Material;
                }
                renderer->Draw(material, rendering_object);
            }
        }
        if (node->Height > 0) {
            for (int i = 0; i < 3; ++i) {
                if (!node->Children[i]) {
                    break;
                }
                Render(group, node, renderer, classification);
            }
        }
    }

    void RenderingTree::ClearState(FeatureInfo* feature_info, int clear_state) {
        for (int i = 0; i < feature_info->Fragments.GetCount(); ++i) {
            feature_info->Fragments.Get(i)->SetState(clear_state);
        }
        if (feature_info->LeftChild) {
            ClearState(feature_info->LeftChild, clear_state);
        }
        if (feature_info->RightChild) {
            ClearState(feature_info->RightChild, clear_state);
        }
    }

    int RenderingTree::Compare(FeatureInfo* feature_info, Feature* feature, const Array<Feature*>& path) {
        if (feature_info->Feature < feature) {
            return -1;
        }
        if (feature_info->Feature > feature) {
            return 1;
        }
        int c1 = feature_info->Path.GetCount();
        int c2 = path.GetCount();
        int c = c1 < c2 ? c1 : c2;
        for (int i = 0; i < c; ++i) {
            Feature* reference1 = feature_info->Path.Get(i);
            Feature* reference2 = path.Get(i);
            if (reference1 < reference2) {
                return -1;
            }
            if (reference1 > reference2) {
                return 1;
            }
        }
        if (c1 < c2) {
            return -1;
        }
        if (c1 > c2) {
            return 1;
        }
        return 0;
    }

    int RenderingTree::Compare(RenderingGroup* rendering_group, RenderingMaterial* material, bool is_classification_enabled) {
        int n = rendering_group->Material->Compare(material);
        if (n != 0) {
            return n;
        }
        if (rendering_group->IsClassificationEnabled < is_classification_enabled) {
            return -1;
        }
        if (rendering_group->IsClassificationEnabled > is_classification_enabled) {
            return 1;
        }
        return 0;
    }

    void RenderingTree::SetFeatureInfoDirty(Feature* feature) {
        if (!m_feature_info_root) {
            return;
        }
        FeatureInfo* feature_info = FindFeatureInfoRoot(m_feature_info_root, feature);
        if (!feature_info) {
            return;
        }
        SetFeatureInfoDirty(feature_info, feature);
        if (m_dirty_level < 1) {
            m_dirty_level = 1;
        }
    }

    FeatureInfo* RenderingTree::FindFeatureInfoRoot(FeatureInfo* feature_info, Feature* feature) {
        if (feature_info->Feature < feature) {
            if (feature_info->LeftChild) {
                return FindFeatureInfoRoot(feature_info->LeftChild, feature);
            }
            return nullptr;
        }
        if (feature_info->Feature > feature) {
            if (feature_info->RightChild) {
                return FindFeatureInfoRoot(feature_info->RightChild, feature);
            }
            return nullptr;
        }
        return feature_info;
    }

    void RenderingTree::SetFeatureInfoDirty(FeatureInfo* feature_info, Feature* feature) {
        if (feature_info->Feature != feature) {
            return;
        }
        if (feature_info->LeftChild) {
            SetFeatureInfoDirty(feature_info->LeftChild, feature);
        }
        if (feature_info->RightChild) {
            SetFeatureInfoDirty(feature_info->RightChild, feature);
        }
        feature_info->IsDirty = true;
        for (int k = 0; k < feature_info->Fragments.GetCount(); ++k) {
            delete feature_info->Fragments.Get(k);
        }
        feature_info->Fragments.Clear();
    }

    FeatureInfo* RenderingTree::GetFeatureInfo(Feature* feature, const Array<Feature*>& path, bool create_not_exist) {
        if (!m_feature_info_root) {
            if (create_not_exist) {
                m_feature_info_root = new FeatureInfo;
                m_feature_info_root->Feature = feature;
                feature->IncRef();
                m_feature_info_root->Path = path;
                for (int i = 0; i < path.GetCount(); ++i) {
                    path.Get(i)->IncRef();
                }
                m_feature_info_root->Order = 0;
                m_feature_info_root->Complexity = 0;
                m_feature_info_root->Group = nullptr;
                m_feature_info_root->IsDirty = true;
                m_feature_info_root->IsRendered = false;
                m_feature_info_root->LeftChild = nullptr;
                m_feature_info_root->RightChild = nullptr;
            }
            return m_feature_info_root;
        }
        return GetFeatureInfo(m_feature_info_root, feature, path, create_not_exist);
    }

    FeatureInfo* RenderingTree::GetFeatureInfo(FeatureInfo* feature_info, Feature* feature, const Array<Feature*>& path, bool create_not_exist) {
        int n = Compare(feature_info, feature, path);
        if (n < 0) {
            FeatureInfo* child_feature_info = feature_info->LeftChild;
            if (!child_feature_info) {
                if (create_not_exist) {
                    child_feature_info = new FeatureInfo;
                    child_feature_info->Feature = feature;
                    feature->IncRef();
                    child_feature_info->Path = path;
                    for (int i = 0; i < path.GetCount(); ++i) {
                        path.Get(i)->IncRef();
                    }
                    child_feature_info->Order = 0;
                    child_feature_info->Complexity = 0;
                    child_feature_info->Group = nullptr;
                    child_feature_info->IsDirty = true;
                    child_feature_info->IsRendered = false;
                    child_feature_info->LeftChild = nullptr;
                    child_feature_info->RightChild = nullptr;
                    feature_info->LeftChild = child_feature_info;
                }
                return child_feature_info;
            }
            return GetFeatureInfo(child_feature_info, feature, path, create_not_exist);
        }
        if (n > 0) {
            FeatureInfo* child_feature_info = feature_info->RightChild;
            if (!child_feature_info) {
                if (create_not_exist) {
                    child_feature_info = new FeatureInfo;
                    child_feature_info->Feature = feature;
                    feature->IncRef();
                    child_feature_info->Path = path;
                    for (int i = 0; i < path.GetCount(); ++i) {
                        path.Get(i)->IncRef();
                    }
                    child_feature_info->Order = 0;
                    child_feature_info->Complexity = 0;
                    child_feature_info->Group = nullptr;
                    child_feature_info->IsDirty = true;
                    child_feature_info->IsRendered = false;
                    child_feature_info->LeftChild = nullptr;
                    child_feature_info->RightChild = nullptr;
                    feature_info->RightChild = child_feature_info;
                }
                return child_feature_info;
            }
            return GetFeatureInfo(child_feature_info, feature, path, create_not_exist);
        }
        return feature_info;
    }

    RenderingGroup* RenderingTree::GetRenderingGroup(RenderingMaterial* material, bool is_classification_enabled, bool create_not_exist) {
        if (!m_rendering_group_root) {
            if (create_not_exist) {
                m_rendering_group_root = new RenderingGroup;
                m_rendering_group_root->Material = material;
                m_rendering_group_root->IsClassificationEnabled = is_classification_enabled;
                m_rendering_group_root->Root = nullptr;
                m_rendering_group_root->LeftChild = nullptr;
                m_rendering_group_root->RightChild = nullptr;
            }
            return m_rendering_group_root;
        }
        return GetRenderingGroup(m_rendering_group_root, material, is_classification_enabled, create_not_exist);
    }

    RenderingGroup* RenderingTree::GetRenderingGroup(RenderingGroup* rendering_group, RenderingMaterial* material, bool is_classification_enabled, bool create_not_exist) {
        int n = Compare(rendering_group, material, is_classification_enabled);
        if (n < 0) {
            RenderingGroup* child_group = rendering_group->LeftChild;
            if (!child_group) {
                if (create_not_exist) {
                    child_group = new RenderingGroup;
                    child_group->Material = material;
                    child_group->IsClassificationEnabled = is_classification_enabled;
                    child_group->Root = nullptr;
                    child_group->LeftChild = nullptr;
                    child_group->RightChild = nullptr;
                    rendering_group->LeftChild = child_group;
                }
                return child_group;
            }
            return GetRenderingGroup(child_group, material, is_classification_enabled, create_not_exist);
        }
        if (n > 0) {
            RenderingGroup* child_group = rendering_group->RightChild;
            if (!child_group) {
                if (create_not_exist) {
                    child_group = new RenderingGroup;
                    child_group->Material = material;
                    child_group->IsClassificationEnabled = is_classification_enabled;
                    child_group->Root = nullptr;
                    child_group->LeftChild = nullptr;
                    child_group->RightChild = nullptr;
                    rendering_group->RightChild = child_group;
                }
                return child_group;
            }
            return GetRenderingGroup(child_group, material, is_classification_enabled, create_not_exist);
        }
        delete material;
        return rendering_group;
    }

    void RenderingTree::FreeFeatureInfos(FeatureInfo* feature_info) {
        if (feature_info->LeftChild) {
            FreeFeatureInfos(feature_info->LeftChild);
        }
        if (feature_info->RightChild) {
            FreeFeatureInfos(feature_info->RightChild);
        }
        for (int i = 0; i < feature_info->Fragments.GetCount(); ++i) {
            delete feature_info->Fragments.Get(i);
        }
        for (int i = 0; i < feature_info->Path.GetCount(); ++i) {
            feature_info->Path.Get(i)->DecRef();
        }
        feature_info->Feature->DecRef();
        delete feature_info;
    }

    void RenderingTree::FreeRenderingGroups(RenderingGroup* rendering_group) {
        if (rendering_group->LeftChild) {
            FreeRenderingGroups(rendering_group->LeftChild);
        }
        if (rendering_group->RightChild) {
            FreeRenderingGroups(rendering_group->RightChild);
        }
        if (rendering_group->Root) {
            FreeRenderingNodes(rendering_group->Root);
        }
        delete rendering_group->Material;
        delete rendering_group;
    }

    void RenderingTree::FreeRenderingNodes(RenderingNode* rendering_node) {
        if (rendering_node->Height > 0) {
            for (int i = 0; i < 3; ++i) {
                RenderingNode* child = (RenderingNode*)rendering_node->Children[i];
                if (!child) {
                    break;
                }
                FreeRenderingNodes(child);
            }
        }
        for (int i = 0; i < rendering_node->RenderingObjects.GetCount(); ++i) {
            rendering_node->RenderingObjects.Get(i)->DecRef();
        }
        delete rendering_node;
    }

}