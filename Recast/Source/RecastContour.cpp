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

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"


static int getCornerHeight(int x, int y, int i, int dir,
						   const rcCompactHeightfield& chf,
						   bool& isBorderVertex)
{
	const rcCompactSpan& s = chf.spans[i];
	int ch = (int)s.y;
	int dirp = (dir+1) & 0x3;
	
	unsigned int regs[4] = {0,0,0,0};
	
    // 注意这里传入来的 xy 并不是顶点。顶点在 xy dir 方向的右边的那个点
    // 以一个顶点为中心的 4 块相邻格子的 span，取这 4 个 span 中高度最高的 span 作为顶点高度
    // 顺便判断了这个顶点是不是区域边界(border)的顶点
    
	// Combine region and area codes in order to prevent
	// border vertices which are in between two areas to be removed.
	regs[0] = chf.spans[i].reg | (chf.areas[i] << 16);
	
	if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
	{
		const int ax = x + rcGetDirOffsetX(dir);
		const int ay = y + rcGetDirOffsetY(dir);
		const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir);
		const rcCompactSpan& as = chf.spans[ai];
		ch = rcMax(ch, (int)as.y);
		regs[1] = chf.spans[ai].reg | (chf.areas[ai] << 16);
		if (rcGetCon(as, dirp) != RC_NOT_CONNECTED)
		{
			const int ax2 = ax + rcGetDirOffsetX(dirp);
			const int ay2 = ay + rcGetDirOffsetY(dirp);
			const int ai2 = (int)chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dirp);
			const rcCompactSpan& as2 = chf.spans[ai2];
			ch = rcMax(ch, (int)as2.y);
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
		}
	}
	if (rcGetCon(s, dirp) != RC_NOT_CONNECTED)
	{
		const int ax = x + rcGetDirOffsetX(dirp);
		const int ay = y + rcGetDirOffsetY(dirp);
		const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dirp);
		const rcCompactSpan& as = chf.spans[ai];
		ch = rcMax(ch, (int)as.y);
		regs[3] = chf.spans[ai].reg | (chf.areas[ai] << 16);
		if (rcGetCon(as, dir) != RC_NOT_CONNECTED)
		{
			const int ax2 = ax + rcGetDirOffsetX(dir);
			const int ay2 = ay + rcGetDirOffsetY(dir);
			const int ai2 = (int)chf.cells[ax2+ay2*chf.width].index + rcGetCon(as, dir);
			const rcCompactSpan& as2 = chf.spans[ai2];
			ch = rcMax(ch, (int)as2.y);
			regs[2] = chf.spans[ai2].reg | (chf.areas[ai2] << 16);
		}
	}

	// Check if the vertex is special edge vertex, these vertices will be removed later.
	for (int j = 0; j < 4; ++j)
	{
		const int a = j;
		const int b = (j+1) & 0x3;
		const int c = (j+2) & 0x3;
		const int d = (j+3) & 0x3;
		
		// The vertex is a border vertex there are two same exterior cells in a row,
		// followed by two interior cells and none of the regions are out of bounds.
        // border 都是沿着轴的，所以这样判断，相当于顶点周围四个 span 全试了一遍
		const bool twoSameExts = (regs[a] & regs[b] & RC_BORDER_REG) != 0 && regs[a] == regs[b];
		const bool twoInts = ((regs[c] | regs[d]) & RC_BORDER_REG) == 0;
		const bool intsSameArea = (regs[c]>>16) == (regs[d]>>16);
		const bool noZeros = regs[a] != 0 && regs[b] != 0 && regs[c] != 0 && regs[d] != 0;
		if (twoSameExts && twoInts && intsSameArea && noZeros)
		{
			isBorderVertex = true;
			break;
		}
	}
	
	return ch;
}

// 获取轮廓顶点集合
static void walkContour(int x, int y, int i,
						rcCompactHeightfield& chf,
						unsigned char* flags, rcIntArray& points)
{
	// Choose the first non-connected edge
	unsigned char dir = 0;
	while ((flags[i] & (1 << dir)) == 0)
		dir++;
	
	unsigned char startDir = dir;
	int starti = i;
	
	const unsigned char area = chf.areas[i];
	
    // 以一个边界 span 作为起始位置，顺时针方向判断它的4条边：
    // 若当前边是区域分界边，则将边的一个顶点加入到轮廓顶点集中，并继续判断下一条边
    // 若当前边不是区域分界边，则移动到与这条边相邻的 span 中（这个 span 是在同一个区域内），重新判断新的 span 的边

    // 关于顶点的判定：span 的某个 dir 方向是区域边界，那么这个边的右边的点就是顶点
    
    // 如果处理 span 超过 40000 还没有处理完区域边，就直接跳出。留待外面 span 循环进行处理
	int iter = 0;
	while (++iter < 40000)
	{
		if (flags[i] & (1 << dir))
		{
			// Choose the edge corner
			bool isBorderVertex = false;
			bool isAreaBorder = false;
			int px = x;
            // 以一个顶点为中心的 4 块相邻格子的span，取这 4 个 span 中高度最高的 span 作为顶点高度
            // 判定顶点是否是区域边界顶点
			int py = getCornerHeight(x, y, i, dir, chf, isBorderVertex);
			int pz = y;
            // 注意这是取的是顶点，span 代表的是格子，case 3 就是自己所以不需要进行处理
			switch(dir)
			{
				case 0: pz++; break;
				case 1: px++; pz++; break;
				case 2: px++; break;
			}
            // r 为邻接区域 id, 是否是边界区域的顶点、是否是区域边界顶点
			int r = 0;
			const rcCompactSpan& s = chf.spans[i];
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
			{
				const int ax = x + rcGetDirOffsetX(dir);
				const int ay = y + rcGetDirOffsetY(dir);
				const int ai = (int)chf.cells[ax+ay*chf.width].index + rcGetCon(s, dir);
				r = (int)chf.spans[ai].reg;
				if (area != chf.areas[ai])
					isAreaBorder = true;
			}
			if (isBorderVertex)
				r |= RC_BORDER_VERTEX;
			if (isAreaBorder)
				r |= RC_AREA_BORDER;
			points.push(px);
			points.push(py);
			points.push(pz);
			points.push(r);
			
			flags[i] &= ~(1 << dir); // Remove visited edges
            // 顺时针旋转
			dir = (dir+1) & 0x3;  // Rotate CW
		}
		else
		{
			int ni = -1;
			const int nx = x + rcGetDirOffsetX(dir);
			const int ny = y + rcGetDirOffsetY(dir);
			const rcCompactSpan& s = chf.spans[i];
			if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
			{
				const rcCompactCell& nc = chf.cells[nx+ny*chf.width];
				ni = (int)nc.index + rcGetCon(s, dir);
			}
			if (ni == -1)
			{
				// Should not happen.
				return;
			}
			x = nx;
			y = ny;
			i = ni;
            // 逆时针旋转回来
			dir = (dir+3) & 0x3;	// Rotate CCW
		}
		
		if (starti == i && startDir == dir)
		{
			break;
		}
	}
}

// (x,z) 点到线段 (px,pz),(qx,qz) 的距离
static float distancePtSeg(const int x, const int z,
						   const int px, const int pz,
						   const int qx, const int qz)
{
    // l点(x,z), p点(px,pz), q点(qx,qz)
	float pqx = (float)(qx - px);
	float pqz = (float)(qz - pz);
	float dx = (float)(x - px);
	float dz = (float)(z - pz);
    // d = |pq|*|pq|
	float d = pqx*pqx + pqz*pqz;
    // 点乘，t = |pq|*|pl|cos(∠qpl)
	float t = pqx*dx + pqz*dz;
	if (d > 0)
		t /= d; // t = l 点落在 pq 上点 k 的比值
    // t < 0:k 超过了起点，直接计算与起点的距离
    // t > 0:k 超过了终点，直接计算与终点的距离
    // 1 > t > 0:k 在线段上，直接计算 k 与 l 的距离
	if (t < 0)
		t = 0;
	else if (t > 1)
		t = 1;
	// px + t*(pqx|pqz) 是 k 点,其实就是 返回的是 k 与 l 的距离平方
	dx = px + t*pqx - x;
	dz = pz + t*pqz - z;
	
	return dx*dx + dz*dz;
}

static void simplifyContour(rcIntArray& points, rcIntArray& simplified,
							const float maxError, const int maxEdgeLen, const int buildFlags)
{
    // 1. 首先不同区域的过渡点需要保留，放到 simplified 中
    // 2. 对于 simplified 中的每个相邻点对(假设记为AB)，检查 points 中位置在AB之间的点，若这些点到 AB 的距离大于某个值 maxError
    // 则将其中距离最远的点加入到 simplified 中，重复这个过程直到所有点距离简化的轮廓线都不超过 maxError
    // 3. maxEdgeLen 需要满足 flag 中有 RC_CONTOUR_TESS_WALL_EDGES|RC_CONTOUR_TESS_AREA_EDGES，也就是仅细分外部边缘或区域之间的边缘。
    // 细化的方式也比较简单，取 simplified 相邻两点，如果起点是区域边界或者是 border 且两点距离超过 maxEdgeLen，取两个点的中间点(point)插入
    // simplified 两点之间
    // 4. simplified 取点的 r 值规则是：边顶点标志取自当前原始点(point)，相邻区域取自下一个点(point)
    //
    // buildFlags 默认值: const int buildFlags = RC_CONTOUR_TESS_WALL_EDGES
    
	// Add initial points.
	bool hasConnections = false;
	for (int i = 0; i < points.size(); i += 4)
	{
		if ((points[i+3] & RC_CONTOUR_REG_MASK) != 0)
		{
			hasConnections = true;
			break;
		}
	}
	
    // simplified 顺序: x,y,z,i(在 points 中的位置)
	if (hasConnections)
	{
        // 不同区域的过渡点需要保留
		// The contour has some portals to other regions.
		// Add a new point to every location where the region changes.
		for (int i = 0, ni = points.size()/4; i < ni; ++i)
		{
			int ii = (i+1) % ni;
            // 判断的点 r 存储的信息:邻接区域 id, 是否是边界区域的顶点、是否是区域边界顶点
			const bool differentRegs = (points[i*4+3] & RC_CONTOUR_REG_MASK) != (points[ii*4+3] & RC_CONTOUR_REG_MASK);
			const bool areaBorders = (points[i*4+3] & RC_AREA_BORDER) != (points[ii*4+3] & RC_AREA_BORDER);
			if (differentRegs || areaBorders)
			{
                // 是邻接区域的交界或者是 border 的交界
				simplified.push(points[i*4+0]);
				simplified.push(points[i*4+1]);
				simplified.push(points[i*4+2]);
				simplified.push(i);
			}
		}
	}
	
	if (simplified.size() == 0)
	{
		// If there is no connections at all,
		// create some initial points for the simplification process.
		// Find lower-left and upper-right vertices of the contour.
        // 如果 simplified size 为 0， 找到顶点中最左下角和最右上角的顶点放入 simplified
		int llx = points[0];
		int lly = points[1];
		int llz = points[2];
		int lli = 0;
		int urx = points[0];
		int ury = points[1];
		int urz = points[2];
		int uri = 0;
		for (int i = 0; i < points.size(); i += 4)
		{
			int x = points[i+0];
			int y = points[i+1];
			int z = points[i+2];
			if (x < llx || (x == llx && z < llz))
			{
				llx = x;
				lly = y;
				llz = z;
				lli = i/4;
			}
			if (x > urx || (x == urx && z > urz))
			{
				urx = x;
				ury = y;
				urz = z;
				uri = i/4;
			}
		}
		simplified.push(llx);
		simplified.push(lly);
		simplified.push(llz);
		simplified.push(lli);
		
		simplified.push(urx);
		simplified.push(ury);
		simplified.push(urz);
		simplified.push(uri);
	}
	
	// Add points until all raw points are within
	// error tolerance to the simplified shape.
    // 添加点，直到所有原始点都在简化形状的误差范围内
	const int pn = points.size()/4;
	for (int i = 0; i < simplified.size()/4; )
	{
		int ii = (i+1) % (simplified.size()/4);
		
		int ax = simplified[i*4+0];
		int az = simplified[i*4+2];
		int ai = simplified[i*4+3];

		int bx = simplified[ii*4+0];
		int bz = simplified[ii*4+2];
		int bi = simplified[ii*4+3];

		// Find maximum deviation from the segment.
        // 找到最大偏差
		float maxd = 0;
		int maxi = -1;
		int ci, cinc, endi;

		// Traverse the segment in lexilogical order so that the
		// max deviation is calculated similarly when traversing
		// opposite segments.
        // 调整遍历段的顺序，以便在遍历相反的段时以类似方式计算最大距离
		if (bx > ax || (bx == ax && bz > az))
		{
			cinc = 1;
			ci = (ai+cinc) % pn;
			endi = bi;
		}
		else
		{
			cinc = pn-1;
			ci = (bi+cinc) % pn;
			endi = ai;
			rcSwap(ax, bx);
			rcSwap(az, bz);
		}
		
		// Tessellate only outer edges or edges between areas.
		if ((points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0 ||
			(points[ci*4+3] & RC_AREA_BORDER))
		{
            // 如果 ci 所在的 point 没有临接区域或者是挨着 border
			while (ci != endi)
			{
                // 计算点与线段的距离
				float d = distancePtSeg(points[ci*4+0], points[ci*4+2], ax, az, bx, bz);
				if (d > maxd)
				{
					maxd = d;
					maxi = ci;
				}
				ci = (ci+cinc) % pn;
			}
		}
		
		
		// If the max deviation is larger than accepted error,
		// add new point, else continue to next segment.
		if (maxi != -1 && maxd > (maxError*maxError))
		{
			// Add space for the new point.
			simplified.resize(simplified.size()+4);
			const int n = simplified.size()/4;
			for (int j = n-1; j > i; --j)
			{
				simplified[j*4+0] = simplified[(j-1)*4+0];
				simplified[j*4+1] = simplified[(j-1)*4+1];
				simplified[j*4+2] = simplified[(j-1)*4+2];
				simplified[j*4+3] = simplified[(j-1)*4+3];
			}
			// Add the point.
			simplified[(i+1)*4+0] = points[maxi*4+0];
			simplified[(i+1)*4+1] = points[maxi*4+1];
			simplified[(i+1)*4+2] = points[maxi*4+2];
			simplified[(i+1)*4+3] = maxi;
		}
		else
		{
			++i;
		}
	}
	
	// Split too long edges.
    // buildFlags 传入为默认值 RC_CONTOUR_TESS_WALL_EDGES
	if (maxEdgeLen > 0 && (buildFlags & (RC_CONTOUR_TESS_WALL_EDGES|RC_CONTOUR_TESS_AREA_EDGES)) != 0)
	{
		for (int i = 0; i < simplified.size()/4; )
		{
			const int ii = (i+1) % (simplified.size()/4);
			
			const int ax = simplified[i*4+0];
			const int az = simplified[i*4+2];
			const int ai = simplified[i*4+3];
			
			const int bx = simplified[ii*4+0];
			const int bz = simplified[ii*4+2];
			const int bi = simplified[ii*4+3];
			
			// Find maximum deviation from the segment.
			int maxi = -1;
			int ci = (ai+1) % pn;
			
			// Tessellate only outer edges or edges between areas.
            // 仅细分外部边缘或区域之间的边缘。
			bool tess = false;
			// Wall edges.
			if ((buildFlags & RC_CONTOUR_TESS_WALL_EDGES) && (points[ci*4+3] & RC_CONTOUR_REG_MASK) == 0)
				tess = true;
			// Edges between areas.
			if ((buildFlags & RC_CONTOUR_TESS_AREA_EDGES) && (points[ci*4+3] & RC_AREA_BORDER))
				tess = true;
			
			if (tess)
			{
				int dx = bx - ax;
				int dz = bz - az;
				if (dx*dx + dz*dz > maxEdgeLen*maxEdgeLen)
				{
					// Round based on the segments in lexilogical order so that the
					// max tesselation is consistent regardles in which direction
					// segments are traversed.
                    // 以词汇顺序基于段进行取整，以便最大细分与穿过段的方向一致。
					const int n = bi < ai ? (bi+pn - ai) : (bi - ai);
					if (n > 1)
					{
                        // bx 在 ax 的右方且包含 bx == ax && bz > az
						if (bx > ax || (bx == ax && bz > az))
							maxi = (ai + n/2) % pn;
						else
							maxi = (ai + (n+1)/2) % pn;
					}
				}
			}
			
			// If the max deviation is larger than accepted error,
			// add new point, else continue to next segment.
			if (maxi != -1)
			{
				// Add space for the new point.
				simplified.resize(simplified.size()+4);
				const int n = simplified.size()/4;
				for (int j = n-1; j > i; --j)
				{
					simplified[j*4+0] = simplified[(j-1)*4+0];
					simplified[j*4+1] = simplified[(j-1)*4+1];
					simplified[j*4+2] = simplified[(j-1)*4+2];
					simplified[j*4+3] = simplified[(j-1)*4+3];
				}
				// Add the point.
				simplified[(i+1)*4+0] = points[maxi*4+0];
				simplified[(i+1)*4+1] = points[maxi*4+1];
				simplified[(i+1)*4+2] = points[maxi*4+2];
				simplified[(i+1)*4+3] = maxi;
			}
			else
			{
				++i;
			}
		}
	}
	
	for (int i = 0; i < simplified.size()/4; ++i)
	{
		// The edge vertex flag is take from the current raw point,
		// and the neighbour region is take from the next raw point.
        // 可能是为了平均？
        // 边顶点标志取自当前原始点，相邻区域取自下一个点。
		const int ai = (simplified[i*4+3]+1) % pn;
		const int bi = simplified[i*4+3];
		simplified[i*4+3] = (points[ai*4+3] & (RC_CONTOUR_REG_MASK|RC_AREA_BORDER)) | (points[bi*4+3] & RC_BORDER_VERTEX);
	}
	
}

// 通过计算 verts 中任意相邻两个点与原点组成的三角形的面积和为正负来判定 verts 是顺时针还是逆时针
// 函数返回负值，即代表是逆时针，即代表是 hole。
static int calcAreaOfPolygon2D(const int* verts, const int nverts)
{
	int area = 0;
	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		const int* vi = &verts[i*4];
		const int* vj = &verts[j*4];
        // 计算两个点与原点组成的三角形的面积
        // area 正为顺时针，负为逆时针。
		area += vi[0] * vj[2] - vj[0] * vi[2];
	}
    // 除以二向上取整
	return (area+1) / 2;
}

// TODO: these are the same as in RecastMesh.cpp, consider using the same.
// Last time I checked the if version got compiled using cmov, which was a lot faster than module (with idiv).
inline int prev(int i, int n) { return i-1 >= 0 ? i-1 : n-1; }
inline int next(int i, int n) { return i+1 < n ? i+1 : 0; }

// 矢量 P = ( x1, y1 )，Q = ( x2, y2 )
// 差积 P×Q = x1*y2 - x2*y1
// 若 P × Q > 0 , 则 P 转 Q 是顺时针方向。
// 若 P × Q < 0 , 则 P 转 Q 是逆时针方向。
// 若 P × Q = 0 , 则 P 与 Q 共线，但可能同向也可能反向。
// ab 在 ac 的哪个方向
inline int area2(const int* a, const int* b, const int* c)
{
	return (b[0] - a[0]) * (c[2] - a[2]) - (c[0] - a[0]) * (b[2] - a[2]);
}

//	Exclusive or: true iff exactly one argument is true.
//	The arguments are negated to ensure that they are 0/1
//	values.  Then the bitwise Xor operator may apply.
//	(This idea is due to Michael Baldwin.)
inline bool xorb(bool x, bool y)
{
	return !x ^ !y;
}

// Returns true iff c is strictly to the left of the directed
// line through a to b.
// c 是否在 ab 左边
inline bool left(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) < 0;
}

inline bool leftOn(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) <= 0;
}

inline bool collinear(const int* a, const int* b, const int* c)
{
	return area2(a, b, c) == 0;
}

//	Returns true iff ab properly intersects cd: they share
//	a point interior to both segments.  The properness of the
//	intersection is ensured by using strict leftness.
// 如果 ab 正确地与 cd 相交，则返回真
static bool intersectProp(const int* a, const int* b, const int* c, const int* d)
{
	// Eliminate improper cases.
    // 杜绝不当案例:共线
	if (collinear(a,b,c) || collinear(a,b,d) ||
		collinear(c,d,a) || collinear(c,d,b))
		return false;
	
    // 同时满足 cd 两点在向量 ab 两边和 ab 两点在向量 cd 两边，才说明 ab 与 cd 相交
	return xorb(left(a,b,c), left(a,b,d)) && xorb(left(c,d,a), left(c,d,b));
}

// Returns T iff (a,b,c) are collinear and point c lies
// on the closed segement ab.
// c 在 ab 内部
static bool between(const int* a, const int* b, const int* c)
{
	if (!collinear(a, b, c))
		return false;
	// If ab not vertical, check betweenness on x; else on y.
	if (a[0] != b[0])
		return	((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]));
	else
		return	((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]));
}

// Returns true iff segments ab and cd intersect, properly or improperly.
static bool intersect(const int* a, const int* b, const int* c, const int* d)
{
	if (intersectProp(a, b, c, d))
		return true;
	else if (between(a, b, c) || between(a, b, d) ||
			 between(c, d, a) || between(c, d, b))
		return true;
	else
		return false;
}

static bool vequal(const int* a, const int* b)
{
	return a[0] == b[0] && a[2] == b[2];
}

// 判断线段 d0d1 是否与 verts 中的相邻点组成的线段相交，排除 i 附近点
static bool intersectSegCountour(const int* d0, const int* d1, int i, int n, const int* verts)
{
	// For each edge (k,k+1) of P
	for (int k = 0; k < n; k++)
	{
		int k1 = next(k, n);
		// Skip edges incident to i.
		if (i == k || i == k1)
			continue;
		const int* p0 = &verts[k * 4];
		const int* p1 = &verts[k1 * 4];
		if (vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1))
			continue;
		
		if (intersect(d0, d1, p0, p1))
			return true;
	}
	return false;
}

// 如果对角线（pj，i）在i端点附近严格位于多边形内部，则返回true。
static bool	inCone(int i, int n, const int* verts, const int* pj)
{
    // 判断是否在锥形内
    // 两种情况
    // 1. i-1,i 逆时针转向 i-1,i+1
    // 2. i-1,i 顺时针转向 i-1,i+1
    // 这两种情况下如果 j 点在 ∠(i-1)(i)(i+1)顺时针方向所形成的夹角范围内(不包含与其中一条共线)
    // 返回 true，否则返回 false
	const int* pi = &verts[i * 4];
	const int* pi1 = &verts[next(i, n) * 4];
	const int* pin1 = &verts[prev(i, n) * 4];
	
	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
    // 如果 p[i] 是个凸顶点，i+1 在 i-1 和 i 的左边或者在一条线
	if (leftOn(pin1, pi, pi1))
		return left(pi, pj, pin1) && left(pj, pi, pi1);
    // 假定 (i-1,i,i+1) 不共线, 或者 p[i] 反射？
	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1));
}


static void removeDegenerateSegments(rcIntArray& simplified)
{
    // 删除 xz 平面上相等的相邻顶点，否则三角剖分器会混淆。
	// Remove adjacent vertices which are equal on xz-plane,
	// or else the triangulator will get confused.
	int npts = simplified.size()/4;
	for (int i = 0; i < npts; ++i)
	{
		int ni = next(i, npts);
		
		if (vequal(&simplified[i*4], &simplified[ni*4]))
		{
			// Degenerate segment, remove.
			for (int j = i; j < simplified.size()/4-1; ++j)
			{
				simplified[j*4+0] = simplified[(j+1)*4+0];
				simplified[j*4+1] = simplified[(j+1)*4+1];
				simplified[j*4+2] = simplified[(j+1)*4+2];
				simplified[j*4+3] = simplified[(j+1)*4+3];
			}
			simplified.resize(simplified.size()-4);
			npts--;
		}
	}
}


static bool mergeContours(rcContour& ca, rcContour& cb, int ia, int ib)
{
    // +2 是因为有重复的两个点 ia ib
    // 外轮廓是顺时针，hole 是逆时针
    // 先从外轮廓的起始点 ia 顺时针绕一圈回到 ia，
    // 然后从洞的起始点 ib 逆时针绕一圈回到 ib
	const int maxVerts = ca.nverts + cb.nverts + 2;
	int* verts = (int*)rcAlloc(sizeof(int)*maxVerts*4, RC_ALLOC_PERM);
	if (!verts)
		return false;
	
	int nv = 0;
	
    // 注意下面是小于等于相当于 ia 添加了两次
	// Copy contour A.
	for (int i = 0; i <= ca.nverts; ++i)
	{
		int* dst = &verts[nv*4];
		const int* src = &ca.verts[((ia+i)%ca.nverts)*4];
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		nv++;
	}

    // 注意下面是小于等于相当于 ib 添加了两次
	// Copy contour B
	for (int i = 0; i <= cb.nverts; ++i)
	{
		int* dst = &verts[nv*4];
		const int* src = &cb.verts[((ib+i)%cb.nverts)*4];
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		nv++;
	}
	
	rcFree(ca.verts);
	ca.verts = verts;
	ca.nverts = nv;
	
	rcFree(cb.verts);
	cb.verts = 0;
	cb.nverts = 0;
	
	return true;
}

struct rcContourHole
{
	rcContour* contour;
	int minx, minz, leftmost;
};

struct rcContourRegion
{
	rcContour* outline;
	rcContourHole* holes;
	int nholes;
};

struct rcPotentialDiagonal
{
	int vert;
	int dist;
};

// Finds the lowest leftmost vertex of a contour.
static void findLeftMostVertex(rcContour* contour, int* minx, int* minz, int* leftmost)
{
	*minx = contour->verts[0];
	*minz = contour->verts[2];
	*leftmost = 0;
	for (int i = 1; i < contour->nverts; i++)
	{
		const int x = contour->verts[i*4+0];
		const int z = contour->verts[i*4+2];
		if (x < *minx || (x == *minx && z < *minz))
		{
			*minx = x;
			*minz = z;
			*leftmost = i;
		}
	}
}

static int compareHoles(const void* va, const void* vb)
{
	const rcContourHole* a = (const rcContourHole*)va;
	const rcContourHole* b = (const rcContourHole*)vb;
	if (a->minx == b->minx)
	{
		if (a->minz < b->minz)
			return -1;
		if (a->minz > b->minz)
			return 1;
	}
	else
	{
		if (a->minx < b->minx)
			return -1;
		if (a->minx > b->minx)
			return 1;
	}
	return 0;
}


static int compareDiagDist(const void* va, const void* vb)
{
	const rcPotentialDiagonal* a = (const rcPotentialDiagonal*)va;
	const rcPotentialDiagonal* b = (const rcPotentialDiagonal*)vb;
	if (a->dist < b->dist)
		return -1;
	if (a->dist > b->dist)
		return 1;
	return 0;
}


static void mergeRegionHoles(rcContext* ctx, rcContourRegion& region)
{
    // 参考 https://www.processon.com/view/link/5f717ecde0b34d327963eb72
    
    // 假设有一个 region 里内嵌了两个 region。大 region 的内部就形成了两个洞：计算边界的时候，
    // 大 region 会算出来三个边界，分别是一个 outline，和两个 hole，
    // 注意，这两个 hole 并非两个小 region 的 outline，这是两码事，每个 region 都有自己独立的边界
    
    // 1. 找出每个小 region 的 LeftMost。LeftMost 的含义是最左下的点：先看哪个顶点的 x 最小，当有两个点 x 相同的时候，则会看这两个哪个 z 更小。
	// Sort holes from left to right.
	for (int i = 0; i < region.nholes; i++)
		findLeftMostVertex(region.holes[i].contour, &region.holes[i].minx, &region.holes[i].minz, &region.holes[i].leftmost);
	
	qsort(region.holes, region.nholes, sizeof(rcContourHole), compareHoles);
	
	int maxVerts = region.outline->nverts;
	for (int i = 0; i < region.nholes; i++)
		maxVerts += region.holes[i].contour->nverts;
	
	rcScopedDelete<rcPotentialDiagonal> diags((rcPotentialDiagonal*)rcAlloc(sizeof(rcPotentialDiagonal)*maxVerts, RC_ALLOC_TEMP));
	if (!diags)
	{
		ctx->log(RC_LOG_WARNING, "mergeRegionHoles: Failed to allocated diags %d.", maxVerts);
		return;
	}
	
	rcContour* outline = region.outline;
	
    // 2. 按顺序遍历 hole（也是按照 LeftMost 顺序）遍历大 region 的 outline 上每个点，找到和洞的最佳顶点匹配的点，
    // 匹配规则为：
    //  1. 距离最近越好
    //  2. 但是连线不能和 region 的 outline 上其他边相交
    // 3. 按照两点连线，给 region上 切一刀，合并 outline 和 hole。新的 outline 的起点和路径都发生了变化，注意会有位置相同的不同顶点存在
	// Merge holes into the outline one by one.
	for (int i = 0; i < region.nholes; i++)
	{
		rcContour* hole = region.holes[i].contour;
		
		int index = -1;
		int bestVertex = region.holes[i].leftmost;
		for (int iter = 0; iter < hole->nverts; iter++)
		{
            // 找到潜在的连线。“最佳”顶点必须位于轮廓中3个连续顶点所描述的圆锥体中
			// Find potential diagonals.
			// The 'best' vertex must be in the cone described by 3 cosequtive vertices of the outline.
			// ..o j-1
			//   |
			//   |   * best
			//   |
			// j o-----o j+1
			//         :
			int ndiags = 0;
			const int* corner = &hole->verts[bestVertex*4];
			for (int j = 0; j < outline->nverts; j++)
			{
                // 如果对角线（corner，j）在j端点附近严格位于多边形内部，则返回true。
                // 判断 corner 是否在 outline->verts 顺序的三个点组成的锥形内
                // 注意应该是顺时针 ∠(j-1)(j)(j+1) 从 j-1 到 j+1 所形成的角度
                // 有这个判断是存在可能点 corner 与 j 的连线是穿过轮廓从轮廓外部连接到 j
				if (inCone(j, outline->nverts, outline->verts, corner))
				{
					int dx = outline->verts[j*4+0] - corner[0];
					int dz = outline->verts[j*4+2] - corner[2];
					diags[ndiags].vert = j;
					diags[ndiags].dist = dx*dx + dz*dz;
					ndiags++;
				}
			}
            // 按距离排序的
			// Sort potential diagonals by distance, we want to make the connection as short as possible.
			qsort(diags, ndiags, sizeof(rcPotentialDiagonal), compareDiagDist);
			
			// Find a diagonal that is not intersecting the outline not the remaining holes.
            // 找到一条不与轮廓相交的连线，也不与剩余的洞相交。
            // 排序之后是从小到大，所以找到的第一个就可以跳出
			index = -1;
			for (int j = 0; j < ndiags; j++)
			{
				const int* pt = &outline->verts[diags[j].vert*4];
                // 判断选择出来的 pt 与 corner 这条线段是否与 outline 轮廓相交
                // bug: i=>j
                // https://github.com/recastnavigation/recastnavigation/issues?q=intersectSegCountour
				bool intersect = intersectSegCountour(pt, corner, diags[j].vert, outline->nverts, outline->verts);
                // 判断选择出来的 pt 与 corner 这条线段是否与 holes 轮廓相交
				for (int k = i; k < region.nholes && !intersect; k++)
					intersect |= intersectSegCountour(pt, corner, -1, region.holes[k].contour->nverts, region.holes[k].contour->verts);
				if (!intersect)
				{
                    // 不相交才可以
					index = diags[j].vert;
					break;
				}
			}
			// If found non-intersecting diagonal, stop looking.
			if (index != -1)
				break;
			// All the potential diagonals for the current vertex were intersecting, try next vertex.
			bestVertex = (bestVertex + 1) % hole->nverts;
		}
		
		if (index == -1)
		{
			ctx->log(RC_LOG_WARNING, "mergeHoles: Failed to find merge points for %p and %p.", region.outline, hole);
			continue;
		}
        // 切一刀合并, 注意 outline 上的 index 点与 hole 上的 bestVertex 添加了两次
        // index:outline 上的最佳点，bestVertex: hole 上的最佳点
		if (!mergeContours(*region.outline, *hole, index, bestVertex))
		{
			ctx->log(RC_LOG_WARNING, "mergeHoles: Failed to merge contours %p and %p.", region.outline, hole);
			continue;
		}
	}
}


/// @par
///
/// The raw contours will match the region outlines exactly. The @p maxError and @p maxEdgeLen
/// parameters control how closely the simplified contours will match the raw contours.
///
/// Simplified contours are generated such that the vertices for portals between areas match up.
/// (They are considered mandatory vertices.)
///
/// Setting @p maxEdgeLength to zero will disabled the edge length feature.
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcAllocContourSet, rcCompactHeightfield, rcContourSet, rcConfig
bool rcBuildContours(rcContext* ctx, rcCompactHeightfield& chf,
					 const float maxError, const int maxEdgeLen,
					 rcContourSet& cset, const int buildFlags)
{
	rcAssert(ctx);
	
    // buildFlags 默认值  const int buildFlags = RC_CONTOUR_TESS_WALL_EDGES
    
	const int w = chf.width;
	const int h = chf.height;
	const int borderSize = chf.borderSize;
	
	rcScopedTimer timer(ctx, RC_TIMER_BUILD_CONTOURS);
	
	rcVcopy(cset.bmin, chf.bmin);
	rcVcopy(cset.bmax, chf.bmax);
	if (borderSize > 0)
	{
		// If the heightfield was build with bordersize, remove the offset.
		const float pad = borderSize*chf.cs;
		cset.bmin[0] += pad;
		cset.bmin[2] += pad;
		cset.bmax[0] -= pad;
		cset.bmax[2] -= pad;
	}
	cset.cs = chf.cs;
	cset.ch = chf.ch;
	cset.width = chf.width - chf.borderSize*2;
	cset.height = chf.height - chf.borderSize*2;
	cset.borderSize = chf.borderSize;
	cset.maxError = maxError;
	
	int maxContours = rcMax((int)chf.maxRegions, 8);
	cset.conts = (rcContour*)rcAlloc(sizeof(rcContour)*maxContours, RC_ALLOC_PERM);
	if (!cset.conts)
		return false;
	cset.nconts = 0;
	
	rcScopedDelete<unsigned char> flags((unsigned char*)rcAlloc(sizeof(unsigned char)*chf.spanCount, RC_ALLOC_TEMP));
	if (!flags)
	{
		ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'flags' (%d).", chf.spanCount);
		return false;
	}
	
	ctx->startTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
	
    // flags[i] 中的四位为 1 表示该方向是个边，0 表示该方向有连接, 还有可能表示该 span 没有 region 或者是 border
	// Mark boundaries.
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			const rcCompactCell& c = chf.cells[x+y*w];
			for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
			{
				unsigned char res = 0;
				const rcCompactSpan& s = chf.spans[i];
				if (!chf.spans[i].reg || (chf.spans[i].reg & RC_BORDER_REG))
				{
					flags[i] = 0;
					continue;
				}
				for (int dir = 0; dir < 4; ++dir)
				{
					unsigned short r = 0;
					if (rcGetCon(s, dir) != RC_NOT_CONNECTED)
					{
						const int ax = x + rcGetDirOffsetX(dir);
						const int ay = y + rcGetDirOffsetY(dir);
						const int ai = (int)chf.cells[ax+ay*w].index + rcGetCon(s, dir);
						r = chf.spans[ai].reg;
					}
					if (r == chf.spans[i].reg)
						res |= (1 << dir);
				}
				flags[i] = res ^ 0xf; // Inverse, mark non connected edges.
			}
		}
	}
	
	ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
	
	rcIntArray verts(256);
	rcIntArray simplified(64);
	
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			const rcCompactCell& c = chf.cells[x+y*w];
			for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
			{
				if (flags[i] == 0 || flags[i] == 0xf)
				{
					flags[i] = 0;
					continue;
				}
				const unsigned short reg = chf.spans[i].reg;
				if (!reg || (reg & RC_BORDER_REG))
					continue;
				const unsigned char area = chf.areas[i];
				
				verts.clear();
				simplified.clear();
				
				ctx->startTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
                // 获取 x, y，i span 所在区域的轮廓顶点集 verts
                // 已经处理的 span dir 会去掉对应 flags 中的值
                // 如果处理 span 超过 40000 还没有处理完区域边，就直接跳出。留待外面 span 循环进行处理
				walkContour(x, y, i, chf, flags, verts);
				ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_TRACE);
				
				ctx->startTimer(RC_TIMER_BUILD_CONTOURS_SIMPLIFY);
                // 简化顶点集合
				simplifyContour(verts, simplified, maxError, maxEdgeLen, buildFlags);
                // 删除 xz 平面上相等的相邻顶点，否则三角剖分器会混淆。
				removeDegenerateSegments(simplified);
				ctx->stopTimer(RC_TIMER_BUILD_CONTOURS_SIMPLIFY);
				
				
				// Store region->contour remap info.
				// Create contour.
				if (simplified.size()/4 >= 3)
				{
					if (cset.nconts >= maxContours)
					{
						// Allocate more contours.
						// This happens when a region has holes.
                        // 有洞就多了边界，所以就多了 nconts，需要重新分配
						const int oldMax = maxContours;
						maxContours *= 2;
						rcContour* newConts = (rcContour*)rcAlloc(sizeof(rcContour)*maxContours, RC_ALLOC_PERM);
						for (int j = 0; j < cset.nconts; ++j)
						{
							newConts[j] = cset.conts[j];
							// Reset source pointers to prevent data deletion.
							cset.conts[j].verts = 0;
							cset.conts[j].rverts = 0;
						}
						rcFree(cset.conts);
						cset.conts = newConts;
						
						ctx->log(RC_LOG_WARNING, "rcBuildContours: Expanding max contours from %d to %d.", oldMax, maxContours);
					}
					
					rcContour* cont = &cset.conts[cset.nconts++];
					
                    // 拷贝简化点的数据到 nverts(size) verts(data)
					cont->nverts = simplified.size()/4;
					cont->verts = (int*)rcAlloc(sizeof(int)*cont->nverts*4, RC_ALLOC_PERM);
					if (!cont->verts)
					{
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'verts' (%d).", cont->nverts);
						return false;
					}
					memcpy(cont->verts, &simplified[0], sizeof(int)*cont->nverts*4);
					if (borderSize > 0)
					{
						// If the heightfield was build with bordersize, remove the offset.
						for (int j = 0; j < cont->nverts; ++j)
						{
							int* v = &cont->verts[j*4];
							v[0] -= borderSize;
							v[2] -= borderSize;
						}
					}
					
                    // 拷贝所有边界点的数据到 nrverts(size) rverts(data)
					cont->nrverts = verts.size()/4;
					cont->rverts = (int*)rcAlloc(sizeof(int)*cont->nrverts*4, RC_ALLOC_PERM);
					if (!cont->rverts)
					{
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'rverts' (%d).", cont->nrverts);
						return false;
					}
					memcpy(cont->rverts, &verts[0], sizeof(int)*cont->nrverts*4);
					if (borderSize > 0)
					{
						// If the heightfield was build with bordersize, remove the offset.
						for (int j = 0; j < cont->nrverts; ++j)
						{
							int* v = &cont->rverts[j*4];
							v[0] -= borderSize;
							v[2] -= borderSize;
						}
					}
					
					cont->reg = reg;
					cont->area = area;
				}
			}
		}
	}
	
	// Merge holes if needed.
	if (cset.nconts > 0)
	{
        // 计算所有轮廓顺时针逆时针
		// Calculate winding of all polygons.
		rcScopedDelete<signed char> winding((signed char*)rcAlloc(sizeof(signed char)*cset.nconts, RC_ALLOC_TEMP));
		if (!winding)
		{
			ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'hole' (%d).", cset.nconts);
			return false;
		}
		int nholes = 0;
		for (int i = 0; i < cset.nconts; ++i)
		{
			rcContour& cont = cset.conts[i];
			// If the contour is wound backwards, it is a hole.
            // 如果轮廓向后缠绕(逆时针)，则它是一个孔。通过计算面积的正负值来判定是顺时针还是逆时针，如果是逆时针的话就是个洞
			winding[i] = calcAreaOfPolygon2D(cont.verts, cont.nverts) < 0 ? -1 : 1;
			if (winding[i] < 0)
				nholes++;
		}
		
		if (nholes > 0)
		{
			// Collect outline contour and holes contours per region.
			// We assume that there is one outline and multiple holes.
            // 收集每个区域的轮廓和洞轮廓。
            // 我们假设有一个轮廓和多个洞。
			const int nregions = chf.maxRegions+1;
			rcScopedDelete<rcContourRegion> regions((rcContourRegion*)rcAlloc(sizeof(rcContourRegion)*nregions, RC_ALLOC_TEMP));
			if (!regions)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'regions' (%d).", nregions);
				return false;
			}
			memset(regions, 0, sizeof(rcContourRegion)*nregions);
			
			rcScopedDelete<rcContourHole> holes((rcContourHole*)rcAlloc(sizeof(rcContourHole)*cset.nconts, RC_ALLOC_TEMP));
			if (!holes)
			{
				ctx->log(RC_LOG_ERROR, "rcBuildContours: Out of memory 'holes' (%d).", cset.nconts);
				return false;
			}
			memset(holes, 0, sizeof(rcContourHole)*cset.nconts);
			
			for (int i = 0; i < cset.nconts; ++i)
			{
				rcContour& cont = cset.conts[i];
				// Positively would contours are outlines, negative holes.
                // 正值为区域轮廓，负值为洞。
				if (winding[i] > 0)
				{
					if (regions[cont.reg].outline)
						ctx->log(RC_LOG_ERROR, "rcBuildContours: Multiple outlines for region %d.", cont.reg);
					regions[cont.reg].outline = &cont;
				}
				else
				{
					regions[cont.reg].nholes++;
				}
			}
			int index = 0;
			for (int i = 0; i < nregions; i++)
			{
				if (regions[i].nholes > 0)
				{
					regions[i].holes = &holes[index];
					index += regions[i].nholes;
					regions[i].nholes = 0;
				}
			}
			for (int i = 0; i < cset.nconts; ++i)
			{
				rcContour& cont = cset.conts[i];
				rcContourRegion& reg = regions[cont.reg];
				if (winding[i] < 0)
					reg.holes[reg.nholes++].contour = &cont;
			}
			
			// Finally merge each regions holes into the outline.
            // 组装好数据 merge region holes
			for (int i = 0; i < nregions; i++)
			{
				rcContourRegion& reg = regions[i];
				if (!reg.nholes) continue;
				
				if (reg.outline)
				{
                    // 对于每个洞，取洞上的一个点和外面轮廓一个点，沿着两点将洞的轮廓添加进现有轮廓中，组成一个更大的轮廓
                    // 其实就是消除洞
					mergeRegionHoles(ctx, reg);
				}
				else
				{
					// The region does not have an outline.
					// This can happen if the contour becaomes selfoverlapping because of
					// too aggressive simplification settings.
                    // 该地区没有轮廓。如果轮廓由于过于激进的简化设置而变得自重叠，则可能会发生这种情况。
					ctx->log(RC_LOG_ERROR, "rcBuildContours: Bad outline for region %d, contour simplification is likely too aggressive.", i);
				}
			}
		}
		
	}
	
	return true;
}
