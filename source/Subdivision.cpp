/****************************************************************************
Copyright 2017, Cognitics Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy 
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.
****************************************************************************/
#include "Subdivision.h"

namespace ctl {

//!    Allocate more memory for Vertices
	void Subdivision::allocVerts(void)
	{
		int start = int(verts.size());
		Block<Vertex>* _verts = new Block<Vertex>(resizeIncriment);
		verts.PushBlock(_verts);
		vertFlags.resize(start + resizeIncriment, false);
		
		for (unsigned int i = 0; i < resizeIncriment; i++)
			_verts->data[i].id = start + i;
	}

//!    Allocate more memory for Edges
	void Subdivision::allocEdges(void)
	{
		int start = getMaxEdges();
		Block<Edge>* _edges = new Block<Edge>(4*resizeIncriment);
		edges.PushBlock(_edges);
		edgeFlags.resize(start + resizeIncriment, false);

		for (size_t i = 0; i < resizeIncriment; i++)
		{
			_edges->data[ 4*i + 0 ].id = static_cast<ID>(start + i);
			_edges->data[ 4*i + 1 ].id = static_cast<ID>(start + i);
			_edges->data[ 4*i + 2 ].id = static_cast<ID>(start + i);
			_edges->data[ 4*i + 3 ].id = static_cast<ID>(start + i);
			_edges->data[ 4*i + 0 ].num = 0;
			_edges->data[ 4*i + 1 ].num = 1;
			_edges->data[ 4*i + 2 ].num = 2;
			_edges->data[ 4*i + 3 ].num = 3;
		}
	}

	Subdivision::Subdivision(int _resizeIncriment)
	{
		resizeIncriment = _resizeIncriment;
		allocVerts();
		allocEdges();
	}

	int Subdivision::getMaxVerts(void) const
	{
		return int(vertFlags.size());
	}

	int Subdivision::getNumVerts(void) const {
		return vertIDGen.numUsedIDs();
	}

	Vertex* Subdivision::getVertex(ID id) {
		return vertFlags[id] ? verts[id] : NULL;
	}

	Vertex* Subdivision::getRandomVertex(void)
	{
		int max = getNumVerts();
		if (max == 0) return NULL;

		int n = rand()%max;
		for (int i = 0; i < max; i++)
		{
			Vertex* vert = getVertex(n);
			if (vert)
				return vert;
			else
				n = (n+1)%max;
		}
		return NULL;
	}

	Vertex* Subdivision::CreateVertex(Point point)
	{
		ID id = vertIDGen.getID();
		if (int(id) >= getMaxVerts())
			allocVerts();
		vertFlags[id] = true;
		Vertex* vert = verts[id];
		vert->clearEdge();
		vert->point = point;
		return vert;
	}

	void Subdivision::RemoveVertex(Vertex* vert)
	{
		if (!vert) return;

		vertFlags[vert->getID()] = false;
		vertIDGen.freeID(vert->getID());

		Edge* edge = vert->getEdges();
		if (edge)
		{
			Edge* next = edge;
			do
			{
				next->setOrg(NULL);
				next = next->Onext();
			} while (next != edge);
		}
	}

	int Subdivision::getMaxEdges(void) const {
		return int(edgeFlags.size());
	}

	int Subdivision::getNumEdges(void) const {
		return edgeIDGen.numUsedIDs();
	}

	Edge* Subdivision::getEdge(ID id) {
		return edgeFlags[id] ? edges[4*id] : NULL;
	}

	Edge* Subdivision::getRandomEdge(void)
	{
		int max = getNumEdges();
		if (max == 0) return NULL;

		int n = rand()%max;
		for (int i = 0; i < max; i++)
		{
			Edge* edge = getEdge(n);
			if (edge)
				return edge;
			else
				n = (n+1)%max;
		}
		return NULL;
	}

	Edge* Subdivision::CreateEdge(Vertex* a, Vertex* b)
	{
		ID id = edgeIDGen.getID();
		if (int(id) >= getMaxEdges())
			allocEdges();
		edgeFlags[id] = true;
		Edge* e0 = edges[ 4*id + 0 ];
		Edge* e1 = edges[ 4*id + 1 ];
		Edge* e2 = edges[ 4*id + 2 ];
		Edge* e3 = edges[ 4*id + 3 ];
		e0->next = e0;
		e1->next = e3;
		e2->next = e2;
		e3->next = e1;
		e0->setEndPoints(a,b);
		e1->point = NULL;
		e3->point = NULL;
		return e0;
	}

	void Subdivision::RemoveEdge(Edge* edge)
	{
		if (!edge) return;

		if (edge->Org()) edge->Org()->removeEdge(edge);
		if (edge->Dest()) edge->Dest()->removeEdge(edge->Sym());
		Edge::Splice(edge,edge->Oprev());
		Edge::Splice(edge->Sym(),edge->Sym()->Oprev());

		edgeIDGen.freeID( edge->getID() );
		edgeFlags[edge->getID()] = false;
	}

	Edge* Subdivision::Connect(Edge* a, Edge* b)
	{
		Edge* edge = CreateEdge( a->Dest() , b->Org() );
		if (!edge) return NULL;

		Splice(edge,a->Lnext());
		Splice(edge->Sym(),b);

		return edge;
	}

	void Subdivision::Splice(Edge* a, Edge* b)
	{
		Edge* alpha = a->Onext()->Rot();
		Edge* beta    = b->Onext()->Rot();

		Edge* temp1 = b->Onext();
		Edge* temp2 = a->Onext();
		Edge* temp3 = beta->Onext();
		Edge* temp4 = alpha->Onext();

		a->next        = temp1;
		b->next        = temp2;
		alpha->next = temp3;
		beta->next    = temp4;
	}

}