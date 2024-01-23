#pragma once
#include "triMesh.h"


namespace TRIANGLE_MESH
{ 
	template <typename vertType>
	struct versVecStatic 
	{
		std::vector<vertType>& vertices;
		size_t versCountValid;
		size_t versSize;
		versVecStatic() = delete;
		versVecStatic(std::vector<vertType>& vers0) : vertices(vers0), versCountValid(vers0.size()), versSize(vers0.size())
		{}

		void upDate()
		{
			this->vertices.resize(versCountValid); 
			versSize = versCountValid; 
		}

		std::vector<vertType>& getData() 
		{
			return this->vertices;
		}

		const std::vector<vertType>& getData() const 
		{
			return this->vertices;
		}

		bool copyData(const std::vector<vertType>& vec0)
		{
			size_t versCount = vec0.size();
			if (this->versSize < versCount)
				return false;
			this->versCountValid = versCount;
			for (size_t i = 0; i < versCount; ++i)
				this->vertices[i] = vec0[i];
			return true;
		}
	};

	////////////////////////////////////////////////////////////////////////////////////////////// triMeshStatic静态三角网格类：
	template <typename TV, typename TI>
	struct triMeshStatic							// 构造时即分配内存，此后其数据所占用的内存空间只会减少不会增加；
	{
		triMesh<TV, TI>& tMesh;
		size_t versCountValid;					// 当前实际的顶点数；
		size_t trisCountValid;
		size_t versSize;								// 顶点数据分配的内存空间可存储的顶点数；
		size_t trisSize;

		triMeshStatic() = delete;
		triMeshStatic(triMesh<TV, TI>& tMesh0) : tMesh(tMesh0), versCountValid(tMesh0.vertices.size()),\
			trisCountValid(tMesh0.triangles.size()), versSize(tMesh0.vertices.size()), trisSize(tMesh0.triangles.size())
		{}

		void upDate()
		{ 
			this->tMesh.vertices.resize(versCountValid);
			this->tMesh.triangles.resize(trisCountValid);
			versSize = versCountValid;
			trisSize = trisCountValid; 
		}

		void clear() 
		{
			this->tMesh.vertices.clear();
			this->tMesh.triangles.clear();
			this->versCountValid = 0;
			this->trisCountValid = 0;
			this->versSize = 0;
			this->trisSize = 0;
		}		 

		triMesh<TV, TI>& getData() 
		{
			return this->tMesh;
		}

		const triMesh<TV, TI>& getData() const 
		{
			return this->tMesh;
		}

		bool copyData(const triMesh<TV, TI>& mesh0) 
		{
			size_t versCount = mesh0.vertices.size();
			size_t trisCount = mesh0.triangles.size();
			if (this->versSize < versCount || this->trisSize < trisCount)
				return false;
			this->versCountValid = versCount; 
			this->trisCountValid = trisCount;
			for (size_t i = 0; i < versCount; ++i)
				this->tMesh.vertices[i] = mesh0.vertices[i];
			for (size_t i = 0; i < trisCount; ++i)
				this->tMesh.triangles[i] = mesh0.triangles[i];
			return true;
		}
	};	
}
using namespace TRIANGLE_MESH;


template <typename TV, typename TI>
triMeshStatic<TV, TI>& getStaticMesh(triMesh<TV, TI>& tMesh0) 
{
	return triMeshStatic<TV, TI>(tMesh0);
}


using triMeshStaticF = triMeshStatic<float, int>;
using triMeshStaticD =triMeshStatic<double, int>;
using versVecStaticF = versVecStatic<verF>;
using versVecStaticD = versVecStatic<verD>;