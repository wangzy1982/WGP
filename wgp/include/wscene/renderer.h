/*
    Original Author: Zuoyuan Wang
    Copyright (c) 2024 Zuoyuan Wang
*/
#ifndef _WGP_SCENE_RENDERER_
#define _WGP_SCENE_RENDERER_

#include "wbase.h"
#include "drawing.h"
#include "wstd/matrix.h"

namespace wgp {

    class RenderingMaterial;
    class RenderingObjectFragment;

    class WGP_API RenderingObject : public RefObject {
    public:
        RenderingObject(int classification, RenderingMaterial* material);
        RenderingMaterial* GetMaterial() const;
        int GetClassification() const;
        virtual ~RenderingObject();
        virtual RenderingObjectFragment* NewFragment() = 0;
        virtual RenderingObjectFragment* Merge(RenderingObject* render_object) = 0;
        virtual RenderingObjectFragment* Merge(RenderingObjectFragment* fragment) = 0;
        virtual void Erase(RenderingObjectFragment* fragment) = 0;
    protected:
        RenderingMaterial* m_material;
        int m_classification;
    };

    class WGP_API RenderingObjectFragment {
    public:
        RenderingObjectFragment(RenderingObject* rendering_object);
        virtual ~RenderingObjectFragment();
        RenderingObject* GetRenderingObject() const;
        virtual int GetClassification() const;
        virtual RenderingObject* NewRenderingObject() const = 0;
    protected:
        RenderingObject* m_rendering_object;
    };

    class WGP_API NullRenderingObjectFragment : public RenderingObjectFragment {
    public:
        NullRenderingObjectFragment(int classification);
        virtual int GetClassification() const;
        virtual RenderingObject* NewRenderingObject() const;
    protected:
        int m_classification;
    };

    class WGP_API RenderingMaterial {
    public:
        virtual ~RenderingMaterial() {}
        virtual int Compare(RenderingMaterial* material) = 0;
        virtual void Render(RenderingObject* rendering_object) = 0;
    };

    struct WGP_API RenderingNode {
        int Height;
        Interval3d Box;
        int Complexity;
        void* Children[3];
        bool IsDirty;
        Array<RenderingObject*> RenderingObjects;
    };
    
    struct WGP_API RenderingGroup {
        RenderingMaterial* Material;
        bool IsClassificationEnabled;
        RenderingNode* Root;
        Array<int> Classifications;
        RenderingGroup* LeftChild;
        RenderingGroup* RightChild;
    };

    struct WGP_API FeatureInfo {
        Feature* Feature;
        Array<wgp::Feature*> Path;
        int Order;
        Interval3d Box;
        int Complexity;
        RenderingGroup* Group;
        bool IsDirty;
        bool IsRendered;
        Array<RenderingObjectFragment*> Fragments;
        FeatureInfo* LeftChild;
        FeatureInfo* RightChild;
    };

    class WGP_API RenderingTree : public DrawingObserver {
    public:
        RenderingTree(Model* model, bool is_order_affected_rendering, int complexity);
        virtual ~RenderingTree();
        virtual void Notify(const Array<CommandLog*>& logs);
        void GetRenderingObjects(int classification, Array<RenderingObject*>& rendering_objects);
        void RemoveRenderingObjects(int classification);
    protected:
        virtual int GetDirtyLevel(CommandLog* log) = 0;
        virtual bool IsRenderingFeature(Feature* feature) = 0;
        virtual void SortFeatures(Array<Feature*>& features) = 0;
        virtual Model* GetReferenceModel(Feature* feature, Matrix4x4& matrix) = 0;
        virtual void Calculate(FeatureInfo* feature_info, RenderingMaterial*& material, bool& is_classification_enabled, Interval3d& box, int& complexity) = 0;
        virtual void BuildRenderingObject(FeatureInfo* feature_info, int classification, Array<RenderingObject*>& rendering_objects) = 0;
    protected:
        void Update();
        void UpdateDirty1(FeatureInfo* feature_info);
        void UpdateDirty2(Array<Feature*>& path, const Matrix4x4& matrix, int& order, Model* model);
        void RemoveUnrenderedFeatureInfo(FeatureInfo** feature_info_address);
        void UpdateDirtyFeatureInfo(FeatureInfo* feature_info);
        void RemoveFeatureInfo(FeatureInfo* feature_info);
        void AddFeatureInfoToGroup(FeatureInfo* feature_info);
        RenderingNode* AddFeatureInfoToGroup(RenderingNode* node, FeatureInfo* feature_info);
        RenderingNode* ReinsertRenderingNodeToGroup(RenderingNode* node, RenderingNode* reinsert_node);
        bool RemoveFeatureInfoFromGroup(FeatureInfo* feature_info);
        bool RemoveFeatureInfoFromGroup(RenderingNode* node, FeatureInfo* feature_info, 
            Array<RenderingNode*>& reinsert_nodes, FeatureInfo*& reinsert_feature_info);
        void BuildRenderingObjects(RenderingGroup* group);
        void BuildRenderingObjects(RenderingGroup* group, RenderingNode* node);
        void BuildRenderingObjects(RenderingGroup* group, RenderingNode* node, Array<RenderingObject*>& rendering_objects);
        void BuildRenderingObjects(int classification, FeatureInfo* feature_info, Array<RenderingObject*>& rendering_objects);
        void AppendRenderingObjects(RenderingGroup* group, int classification);
        void AppendRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node);
        void AppendRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node, Array<RenderingObject*>& rendering_objects);
        void AppendRenderingObjects(int classification, FeatureInfo* feature_info, Array<RenderingObject*>& rendering_objects);
        void RemoveRenderingObjects(RenderingGroup* group, int classification);
        void RemoveRenderingObjects(RenderingGroup* group, int classification, RenderingNode* node);
        void RemoveRenderingObjects(int classification, FeatureInfo* feature_info);
        int Compare(FeatureInfo* feature_info, Feature* feature, const Array<Feature*>& path);
        int Compare(RenderingGroup* rendering_group, RenderingMaterial* material, bool is_classification_enabled);
        void SetFeatureInfoDirty(Feature* feature);
        FeatureInfo* FindFeatureInfoRoot(FeatureInfo* feature_info, Feature* feature);
        void SetFeatureInfoDirty(FeatureInfo* feature_info, Feature* feature);
        FeatureInfo* GetFeatureInfo(Feature* feature, const Array<Feature*>& path, bool create_not_exist);
        FeatureInfo* GetFeatureInfo(FeatureInfo* feature_info, Feature* feature, const Array<Feature*>& path, bool create_not_exist);
        RenderingGroup* GetRenderingGroup(RenderingMaterial* material, bool is_classification_enabled, bool create_not_exist);
        RenderingGroup* GetRenderingGroup(RenderingGroup* rendering_group, RenderingMaterial* material, bool is_classification_enabled, bool create_not_exist);
        void FreeFeatureInfos(FeatureInfo* feature_info);
        void FreeRenderingGroups(RenderingGroup* rendering_group);
        void FreeRenderingNodes(RenderingNode* rendering_node);
    private:
        Model* m_model;
        bool m_is_order_affected_rendering;
        FeatureInfo* m_feature_info_root;
        RenderingGroup* m_rendering_group_root;
        int m_complexity;
        int m_dirty_level;  //0: Not dirty;   1: Changed;    2: Added or removed
    };

}

#endif