//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef CHUNKYTRIMESH_H
#define CHUNKYTRIMESH_H

struct rcChunkyTriMeshNode
{
	float bmin[2]; //包含区域内所有三角形的包围盒
	float bmax[2];
	int i;  //i为正表示是叶子节点；i为负表示是非叶子节点，大小等于在节点数组内的覆盖范围
	int n; //节点内包含的三角形个数
};

struct rcChunkyTriMesh
{
	inline rcChunkyTriMesh() : nodes(0), nnodes(0), tris(0), ntris(0), maxTrisPerChunk(0) {};
	inline ~rcChunkyTriMesh() { delete [] nodes; delete [] tris; }

	rcChunkyTriMeshNode* nodes; // 树的节点
	int nnodes; // 节点的数量
	int* tris;  // 三角形数组 i*3 i[0-2] verts[i[0-2]]
	int ntris;  // 三角形的数量
	int maxTrisPerChunk; // 叶子节点中三角形最大数量

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	rcChunkyTriMesh(const rcChunkyTriMesh&);
	rcChunkyTriMesh& operator=(const rcChunkyTriMesh&);
};

/// Creates partitioned triangle mesh (AABB tree),
/// where each node contains at max trisPerChunk triangles.
bool rcCreateChunkyTriMesh(const float* verts, const int* tris, int ntris,
						   int trisPerChunk, rcChunkyTriMesh* cm);

/// Returns the chunk indices which overlap the input rectable.
int rcGetChunksOverlappingRect(const rcChunkyTriMesh* cm, float bmin[2], float bmax[2], int* ids, const int maxIds);

/// Returns the chunk indices which overlap the input segment.
int rcGetChunksOverlappingSegment(const rcChunkyTriMesh* cm, float p[2], float q[2], int* ids, const int maxIds);


#endif // CHUNKYTRIMESH_H
