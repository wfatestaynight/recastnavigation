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
#include <stdio.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"

inline bool overlapBounds(const float* amin, const float* amax, const float* bmin, const float* bmax)
{
	bool overlap = true;
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
	return overlap;
}

inline bool overlapInterval(unsigned short amin, unsigned short amax,
							unsigned short bmin, unsigned short bmax)
{
	if (amax < bmin) return false;
	if (amin > bmax) return false;
	return true;
}


static rcSpan* allocSpan(rcHeightfield& hf)
{
	// If running out of memory, allocate new page and update the freelist.
	if (!hf.freelist || !hf.freelist->next)
	{
		// Create new page.
		// Allocate memory for the new pool.
		rcSpanPool* pool = (rcSpanPool*)rcAlloc(sizeof(rcSpanPool), RC_ALLOC_PERM);
		if (!pool) return 0;

		// Add the pool into the list of pools.
		pool->next = hf.pools;
		hf.pools = pool;
		// Add new items to the free list.
		rcSpan* freelist = hf.freelist;
		rcSpan* head = &pool->items[0];
		rcSpan* it = &pool->items[RC_SPANS_PER_POOL];
		do
		{
			--it;
			it->next = freelist;
			freelist = it;
		}
		while (it != head);
		hf.freelist = it;
	}
	
	// Pop item from in front of the free list.
	rcSpan* it = hf.freelist;
	hf.freelist = hf.freelist->next;
	return it;
}

static void freeSpan(rcHeightfield& hf, rcSpan* ptr)
{
	if (!ptr) return;
	// Add the node in front of the free list.
	ptr->next = hf.freelist;
	hf.freelist = ptr;
}

static bool addSpan(rcHeightfield& hf, const int x, const int y,
					const unsigned short smin, const unsigned short smax,
					const unsigned char area, const int flagMergeThr)
{
	
    // flagMergeThr 为参数 walkableClimb Max Climb 最大攀爬高度
	int idx = x + y*hf.width;
	
	rcSpan* s = allocSpan(hf);
	if (!s)
		return false;
	s->smin = smin;
	s->smax = smax;
	s->area = area;
	s->next = 0;
	
	// Empty cell, add the first span.
	if (!hf.spans[idx])
	{
		hf.spans[idx] = s;
		return true;
	}
	rcSpan* prev = 0;
	rcSpan* cur = hf.spans[idx];
	
	// Insert and merge spans.
	while (cur)
	{
		if (cur->smin > s->smax)
		{
			// Current span is further than the new span, break.
			break;
		}
		else if (cur->smax < s->smin)
		{
			// Current span is before the new span advance.
			prev = cur;
			cur = cur->next;
		}
		else
		{
            // cur 与 s 有重合的部分
			// Merge spans.
			if (cur->smin < s->smin)
				s->smin = cur->smin;
			if (cur->smax > s->smax)
				s->smax = cur->smax;
			
            // smax 是最高处的差在可攀爬范围内就合并 flag
            // 如果同一个 span，有的面是可行走，有的是不可行走，在合并之后会变成可行走
			// Merge flags.
			if (rcAbs((int)s->smax - (int)cur->smax) <= flagMergeThr)
				s->area = rcMax(s->area, cur->area);
			
			// Remove current span.
			rcSpan* next = cur->next;
			freeSpan(hf, cur);
			if (prev)
				prev->next = next;
			else
				hf.spans[idx] = next;
			cur = next;
		}
	}
	
    // hf.spans 从链表开始到结束，3d 中的 y 值越来越高
	// Insert new span.
	if (prev)
	{
		s->next = prev->next;
		prev->next = s;
	}
	else
	{
		s->next = hf.spans[idx];
		hf.spans[idx] = s;
	}

	return true;
}

/// @par
///
/// The span addition can be set to favor flags. If the span is merged to
/// another span and the new @p smax is within @p flagMergeThr units
/// from the existing span, the span flags are merged.
///
/// @see rcHeightfield, rcSpan.
bool rcAddSpan(rcContext* ctx, rcHeightfield& hf, const int x, const int y,
			   const unsigned short smin, const unsigned short smax,
			   const unsigned char area, const int flagMergeThr)
{
	rcAssert(ctx);

	if (!addSpan(hf, x, y, smin, smax, area, flagMergeThr))
	{
		ctx->log(RC_LOG_ERROR, "rcAddSpan: Out of memory.");
		return false;
	}

	return true;
}

// 用一条与某一坐标轴平行的线分割凸多边形成两个凸多边形, 小于 x 的 是 out1
// divides a convex polygons into two convex polygons on both sides of a line
static void dividePoly(const float* in, int nin,
					  float* out1, int* nout1,
					  float* out2, int* nout2,
					  float x, int axis)
{
    // axis 是 0,1,2 对应 xyz, x 为 axis 对应(xyz)的某一个值
	float d[12];
	for (int i = 0; i < nin; ++i)
		d[i] = x - in[i*3+axis];

	int m = 0, n = 0;
    // 循环不重复选择多边形的两个点
	for (int i = 0, j = nin-1; i < nin; j=i, ++i)
	{
		bool ina = d[j] >= 0;
		bool inb = d[i] >= 0;
		if (ina != inb)
		{
            // i 点与 j 点在 x 两边，利用三角形相似计算出 i 与 j 点所连成的线与 x 的交点，并分别放入分割的两个凸多边形点集合中
			float s = d[j] / (d[j] - d[i]);
			out1[m*3+0] = in[j*3+0] + (in[i*3+0] - in[j*3+0])*s;
			out1[m*3+1] = in[j*3+1] + (in[i*3+1] - in[j*3+1])*s;
			out1[m*3+2] = in[j*3+2] + (in[i*3+2] - in[j*3+2])*s;
			rcVcopy(out2 + n*3, out1 + m*3);
			m++;
			n++;
			// add the i'th point to the right polygon. Do NOT add points that are on the dividing line
			// since these were already added above
			if (d[i] > 0)
			{
				rcVcopy(out1 + m*3, in + i*3);
				m++;
			}
			else if (d[i] < 0)
			{
				rcVcopy(out2 + n*3, in + i*3);
				n++;
			}
		}
		else // same side
		{
			// add the i'th point to the right polygon. Addition is done even for points on the dividing line
			if (d[i] >= 0)
			{
				rcVcopy(out1 + m*3, in + i*3);
				m++;
                // 等于 0 就是在分割线上所以要两边都放
				if (d[i] != 0)
					continue;
			}
			rcVcopy(out2 + n*3, in + i*3);
			n++;
		}
	}

	*nout1 = m;
	*nout2 = n;
}



static bool rasterizeTri(const float* v0, const float* v1, const float* v2,
						 const unsigned char area, rcHeightfield& hf,
						 const float* bmin, const float* bmax,
						 const float cs, const float ics, const float ich,
						 const int flagMergeThr)
{
	const int w = hf.width;
	const int h = hf.height;
	float tmin[3], tmax[3];
	const float by = bmax[1] - bmin[1];
	
	// Calculate the bounding box of the triangle.
	rcVcopy(tmin, v0);
	rcVcopy(tmax, v0);
	rcVmin(tmin, v1);
	rcVmin(tmin, v2);
	rcVmax(tmax, v1);
	rcVmax(tmax, v2);
	
    // 前面说过选取的区域不一定是整个 bounding box 的大小， 所以计算的时候跳过这些不需要计算的三角形
	// If the triangle does not touch the bbox of the heightfield, skip the triagle.
	if (!overlapBounds(bmin, bmax, tmin, tmax))
		return true;
	
	// Calculate the footprint of the triangle on the grid's y-axis
	int y0 = (int)((tmin[2] - bmin[2])*ics);
	int y1 = (int)((tmax[2] - bmin[2])*ics);
    // 限制 y0 y1 在 0 到 h-1 之间
    // 这个 y 其实是 3d 中的 z
	y0 = rcClamp(y0, 0, h-1);
	y1 = rcClamp(y1, 0, h-1);
	
	// Clip the triangle into all grid cells it touches.
    // 7:一个三角形在分割时会被切4刀（即一个正方形格子的4条边），所以切完最多可能会变成一个7边形，也就是最多需要存7个顶点
    // 3:对应一个顶点的 x y z 坐标值
    // 4:有4份对应的数据，用来存放输入的图形和分割后的图形
	float buf[7*3*4];
	float *in = buf, *inrow = buf+7*3, *p1 = inrow+7*3, *p2 = p1+7*3;

	rcVcopy(&in[0], v0);
	rcVcopy(&in[1*3], v1);
	rcVcopy(&in[2*3], v2);
	int nvrow, nvIn = 3;
	
	for (int y = y0; y <= y1; ++y)
	{
		// Clip polygon to row. Store the remaining polygon as well
		const float cz = bmin[2] + y*cs;
        // 把 in 表示的多边形以以 z=cz+cs 分割成 inrow(左边) 与 p1(右边)两个多边形
		dividePoly(in, nvIn, inrow, &nvrow, p1, &nvIn, cz+cs, 2);
        // x 轴右边的继续等待下次循环划分
		rcSwap(in, p1);
        // 左边没有数据
		if (nvrow < 3) continue;
		
        // 只算的小于 cz+cs(左边)的数据，左边的是已经划分好的数据在 [cz+cs-1, cz+cs] 之间
		// find the horizontal bounds in the row
		float minX = inrow[0], maxX = inrow[0];
		for (int i=1; i<nvrow; ++i)
		{
			if (minX > inrow[i*3])	minX = inrow[i*3];
			if (maxX < inrow[i*3])	maxX = inrow[i*3];
		}
		int x0 = (int)((minX - bmin[0])*ics);
		int x1 = (int)((maxX - bmin[0])*ics);
		x0 = rcClamp(x0, 0, w-1);
		x1 = rcClamp(x1, 0, w-1);

		int nv, nv2 = nvrow;

		for (int x = x0; x <= x1; ++x)
		{
			// Clip polygon to column. store the remaining polygon as well
			const float cx = bmin[0] + x*cs;
            // 再按照 x 轴划分一次
			dividePoly(inrow, nv2, p1, &nv, p2, &nv2, cx+cs, 0);
            // x 轴右边的继续等待下次循环划分
			rcSwap(inrow, p2);
			if (nv < 3) continue;
			
            // 只算的小于 cx+cs(左边)的数据，左边的是已经划分好的数据在 [cx+cs-1, cx+cs] 之间
			// Calculate min and max of the span.
			float smin = p1[1], smax = p1[1];
			for (int i = 1; i < nv; ++i)
			{
				smin = rcMin(smin, p1[i*3+1]);
				smax = rcMax(smax, p1[i*3+1]);
			}
			smin -= bmin[1];
			smax -= bmin[1];
			// Skip the span if it is outside the heightfield bbox
			if (smax < 0.0f) continue;
			if (smin > by) continue;
			// Clamp the span to the heightfield bbox.
			if (smin < 0.0f) smin = 0;
			if (smax > by) smax = by;
			
			// Snap the span to the heightfield height grid.
			unsigned short ismin = (unsigned short)rcClamp((int)floorf(smin * ich), 0, RC_SPAN_MAX_HEIGHT);
			unsigned short ismax = (unsigned short)rcClamp((int)ceilf(smax * ich), (int)ismin+1, RC_SPAN_MAX_HEIGHT);
			
            // 到这里已经取到了一个格子的 bbox 数据，然后放入高度场中
			if (!addSpan(hf, x, y, ismin, ismax, area, flagMergeThr))
				return false;
		}
	}

	return true;
}

/// @par
///
/// No spans will be added if the triangle does not overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangle(rcContext* ctx, const float* v0, const float* v1, const float* v2,
						 const unsigned char area, rcHeightfield& solid,
						 const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);

	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	if (!rasterizeTri(v0, v1, v2, area, solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
	{
		ctx->log(RC_LOG_ERROR, "rcRasterizeTriangle: Out of memory.");
		return false;
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const int /*nv*/,
						  const int* tris, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	// Rasterize triangles.
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[tris[i*3+0]*3];
		const float* v1 = &verts[tris[i*3+1]*3];
		const float* v2 = &verts[tris[i*3+2]*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const int /*nv*/,
						  const unsigned short* tris, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);

	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	// Rasterize triangles.
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[tris[i*3+0]*3];
		const float* v1 = &verts[tris[i*3+1]*3];
		const float* v2 = &verts[tris[i*3+2]*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

/// @par
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
bool rcRasterizeTriangles(rcContext* ctx, const float* verts, const unsigned char* areas, const int nt,
						  rcHeightfield& solid, const int flagMergeThr)
{
	rcAssert(ctx);
	
	rcScopedTimer timer(ctx, RC_TIMER_RASTERIZE_TRIANGLES);
	
	const float ics = 1.0f/solid.cs;
	const float ich = 1.0f/solid.ch;
	// Rasterize triangles.
	for (int i = 0; i < nt; ++i)
	{
		const float* v0 = &verts[(i*3+0)*3];
		const float* v1 = &verts[(i*3+1)*3];
		const float* v2 = &verts[(i*3+2)*3];
		// Rasterize.
		if (!rasterizeTri(v0, v1, v2, areas[i], solid, solid.bmin, solid.bmax, solid.cs, ics, ich, flagMergeThr))
		{
			ctx->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}
