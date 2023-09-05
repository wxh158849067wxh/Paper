#pragma once
#ifndef PARSENEF_H
#define PARSENEF_H
 // !PARSENEF_H

#include<memory>
#include<vector>
#include<iostream>
#include<assert.h>
#include<unordered_set>
#include<unordered_map>
#include<set>

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL\Nef_3\SNC_iteration.h>
#include<CGAL\circulator.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
namespace wxh {
	namespace internal
	{
		template <class Point, class Vector>
		void newell_single_step_3(const Point& p, const Point& q, Vector& n)
		{
			// Compute normal of the face by using Newell's method: for each edge PQ
			// Nx += (Py - Qy) * (Pz + Qz);
			// Ny += (Pz - Qz) * (Px + Qx);
			// Nz += (Px - Qx) * (Py + Qy);
			n = Vector(n.x() + ((p.y() - q.y())*(p.z() + q.z())),
				n.y() + ((p.z() - q.z())*(p.x() + q.x())),
				n.z() + ((p.x() - q.x())*(p.y() + q.y())));
		}

		template <class Point, class Vector>
		Vector compute_normal_of_face(const std::vector<Point>& points)
		{
			Vector normal(CGAL::NULL_VECTOR);
			unsigned int nb = 0;
			for (std::size_t i = 0; i < points.size(); ++i)
			{
				internal::newell_single_step_3(points[i], points[(i + 1) % points.size()],
					normal);
				++nb;
			}

			assert(nb > 0);
			return (typename CGAL::Kernel_traits<Vector>::Kernel::Construct_scaled_vector_3()
				(normal, 1.0 / nb));
		}

		////////////////////////////////////////////////////////////////
		// Structs to transform any CGAL point/vector into a Local_point/Local_vector
		template<typename K, typename Local_kernel>
		struct Geom_utils
		{
			static typename Local_kernel::Point_3 get_local_point(const typename K::Point_2& p)
			{
				CGAL::Cartesian_converter<K, Local_kernel> converter;
				return typename Local_kernel::Point_3(converter(p.x()), converter(p.y()), 0);
			}
			static typename Local_kernel::Point_3 get_local_point(const typename K::Weighted_point_2& p)
			{
				typename K::Point_2 lp(p);
				return Geom_utils<K, Local_kernel>::get_local_point(lp);
			}
			static typename Local_kernel::Point_3 get_local_point(const typename K::Point_3& p)
			{
				CGAL::Cartesian_converter<K, Local_kernel> converter;
				return converter(p);
			}
			static typename Local_kernel::Point_3 get_local_point(const typename K::Weighted_point_3& p)
			{
				typename K::Point_3 lp(p);
				return Geom_utils<K, Local_kernel>::get_local_point(lp);
			}
			static typename Local_kernel::Vector_3 get_local_vector(const typename K::Vector_2& v)
			{
				CGAL::Cartesian_converter<K, Local_kernel> converter;
				return typename Local_kernel::Vector_3(converter(v.x()), converter(v.y()), 0);
			}
			static typename Local_kernel::Vector_3 get_local_vector(const typename K::Vector_3& v)
			{
				CGAL::Cartesian_converter<K, Local_kernel> converter;
				return converter(v);
			}
			static typename Local_kernel::Ray_2 get_local_ray(const typename K::Ray_2& r)
			{
				CGAL::Cartesian_converter<K, Local_kernel> converter;
				return converter(r);
			}
		};
		template<typename Local_kernel>
		struct Geom_utils<Local_kernel, Local_kernel>
		{
			static typename Local_kernel::Point_3 get_local_point(const typename Local_kernel::Point_2& p)
			{
				return typename Local_kernel::Point_3(p.x(), p.y(), 0);
			}
			static typename Local_kernel::Point_3 get_local_point(const typename Local_kernel::Weighted_point_2& p)
			{
				return typename Local_kernel::Point_3(p.point().x(), p.point().y(), 0);
			}
			static const typename Local_kernel::Point_3 & get_local_point(const typename Local_kernel::Point_3& p)
			{
				return p;
			}
			static typename Local_kernel::Point_3 get_local_point(const typename Local_kernel::Weighted_point_3& p)
			{
				return typename Local_kernel::Point_3(p);
			}
			static typename Local_kernel::Vector_3 get_local_vector(const typename Local_kernel::Vector_2& v)
			{
				return typename Local_kernel::Vector_3(v.x(), v.y(), 0);
			}
			static const typename Local_kernel::Vector_3& get_local_vector(const typename Local_kernel::Vector_3& v)
			{
				return v;
			}
			static const typename Local_kernel::Ray_2& get_local_ray(const typename Local_kernel::Ray_2& r)
			{
				return r;
			}
		};
	}
	template<class Nef_Polyhedron>
	class ParseNef
	{
		typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
		typedef Local_kernel::Point_3  Local_point;
		typedef Local_kernel::Vector_3 Local_vector;
		typedef Local_kernel::Ray_2    Local_ray;
		//typedef typename Container_Vector::value_type Nef_Polyhedron;
		//typedef typename Container_Vector::size_type size_type;
		typedef typename Nef_Polyhedron::Kernel                   Kernel;

		typedef typename Nef_Polyhedron::Halffacet_cycle_const_iterator            Halffacet_cycle_const_iterator;
		typedef typename Nef_Polyhedron::SHalfedge_around_facet_const_circulator   SHalfedge_around_facet_const_circulator;

		typedef typename Nef_Polyhedron::Shell_entry_const_iterator   Shell_entry_const_iterator;
		typedef typename Nef_Polyhedron::SHalfedge_const_iterator     SHalfedge_const_iterator;
		typedef typename Nef_Polyhedron::Volume_const_iterator        Volume_const_iterator;

		typedef typename Nef_Polyhedron::Vertex_const_handle       Vertex_const_handle;
		typedef typename Nef_Polyhedron::SFace_const_handle        SFace_const_handle;
		typedef typename Nef_Polyhedron::Halfedge_const_handle     Halfedge_const_handle;
		typedef typename Nef_Polyhedron::Halffacet_const_handle    Halffacet_const_handle;
		typedef typename Nef_Polyhedron::SHalfedge_const_handle    SHalfedge_const_handle;
		typedef typename Nef_Polyhedron::SHalfloop_const_handle    SHalfloop_const_handle;
		typedef typename Nef_Polyhedron::Point_3 Point;

	protected:
		// Shortcuts to simplify function calls.
		template<typename KPoint>
		static Local_point get_local_point(const KPoint& p)
		{
			return internal::Geom_utils<typename CGAL::Kernel_traits<KPoint>::Kernel, Local_kernel>::
				get_local_point(p);
		}
		template<typename KVector>
		static Local_vector get_local_vector(const KVector& v)
		{
			return internal::Geom_utils<typename CGAL::Kernel_traits<KVector>::Kernel, Local_kernel>::
				get_local_vector(v);
		}
		Local_vector get_vertex_normal(Vertex_const_handle vh);//计算一个顶点的法向量；
		Local_vector get_face_normal(SHalfedge_const_handle she);//计算一个facet的法向量
	private:
		std::vector<float> m_points;//所有顶点
		std::vector<float> m_segments;//所有的线段
		std::vector<float> m_trigulations;//进行了三角剖分的顶点
		std::vector<Local_point> m_faces;//一个facet的所有顶点，未进行三角剖分
		std::vector<std::size_t> m_indexs;//为了drawelemet
		std::vector<float>m_normal_vertex;//每个顶点的法向量
		std::vector<Local_point> m_all_points;//所有的顶点
		std::unordered_map<Local_point, std::size_t> m_map;//为了记录每一个顶点的index，为了得到m_indexs;
		CGAL::Bbox_3 m_bbox;
		std::size_t index = 0;
	public:
		ParseNef() = delete;
		ParseNef(Nef_Polyhedron& cv);
		ParseNef(std::vector<Nef_Polyhedron>& cv);
		~ParseNef() = default;
		void add_point(const Point& p1);
		void add_points(const Local_point& p1);
		void add_segment(const Point& p1, const Point& p2);
		void add_point_in_triangulation(const Local_point& p1);
		void add_point_in_face(const Point& p1);
		//void compute_componet();//解析Nef多面体
		void triangular_face_end_internal(const Local_vector& normal);//当一个facet本身为三角形(m_faces的顶点等于3)不需要进行三角剖分
		bool is_facet_convex(const std::vector<Local_point>& facet, const Local_vector& normal);//判断一个facet是否为凸多边形
		void triangulate(Local_vector& normal, std::vector<Local_point>& vc);
		void convex_quadrangular_face_end_internal(const Local_vector& normal);//当一个facet只有四条边(即m_faces只有四个顶点时)进行三角剖分
		void convex_face_end_internal(const Local_vector& normal);//当一个facet为凸且边数大于4(即m_faces的顶点大于4时)进行三角剖分
		void nonconvex_face_end_internal(const Local_vector& normal);//对凹多边形进行三角剖分
	protected:
		class Nef_Visitor
		{
		public:
			Nef_Visitor(ParseNef<Nef_Polyhedron> *v)
				: n_faces(0), n_edges(0), viewer(v) {}

			void visit(Vertex_const_handle vh) {

				if (vertex_set.count(vh) == 0)
				{
					vertex_set.insert(vh);
					viewer->add_point(vh->point());
				}

			}

			void visit(Halffacet_const_handle facet)
			{
				Halffacet_const_handle f = facet->twin();
				if (face_set.count(f) == 1 || face_set.count(facet) == 1)
				{
					return;
				}
				//if (facets_done.find(f) != facets_done.end() ||
				//	facets_done.find(opposite_facet) != facets_done.end()) {
				//	return;
				//}

				SHalfedge_const_handle se;
				Halffacet_cycle_const_iterator fc;
				fc = f->facet_cycles_begin();

				se = SHalfedge_const_handle(fc); // non-zero if shalfedge is returned
				if (se == 0)
				{ //return if not-shalfedge
					return;
				}
				face_set.insert(facet);
				face_set.insert(f);
				///*CGAL::IO::Color c = viewer.run_color(f);
				//viewer.face_begin(c);*/

				SHalfedge_around_facet_const_circulator hc_start(se);
				SHalfedge_around_facet_const_circulator hc_end(hc_start);
				CGAL_For_all(hc_start, hc_end) {
					Vertex_const_handle vh = hc_start->source()->center_vertex();
					//	//std::cout << vh->point();
					viewer->add_point_in_face(vh->point());
					viewer->get_vertex_normal(vh);
					//}
					////viewer.face_end();
					//facets_done[f] = true;
					//n_faces++;
				}
				Local_vector normal = internal::compute_normal_of_face
					<Local_point, Local_vector>(viewer->m_faces);//得到平面法向量
				viewer->triangulate(normal, viewer->m_faces);

				++n_faces;
				viewer->m_faces.clear();
			}
			void visit(Halfedge_const_handle he)
			{
				Halfedge_const_handle twin = he->twin();
				if (edge_set.find(he) != edge_set.end() ||
					edge_set.find(twin) != edge_set.end())
				{
					// Edge already added
					return;
				}

				viewer->add_segment(he->source()->point(), he->target()->point());
				edge_set.insert(he);
				edge_set.insert(twin);
				n_edges++;
			}

			void visit(SHalfedge_const_handle) {}
			void visit(SHalfloop_const_handle) {}
			void visit(SFace_const_handle) {}
			int n_faces;
			int n_edges;
		protected:
			//std::unordered_map<Halffacet_const_handle, bool> facets_done;
			//std::unordered_map<Halfedge_const_handle, bool> edges_done;
			std::unordered_set<Halffacet_const_handle> face_set;
			std::unordered_set<Halfedge_const_handle> edge_set;
			std::unordered_set<Vertex_const_handle> vertex_set;
			ParseNef<Nef_Polyhedron>* viewer;

		};
	public:
		const CGAL::Bbox_3 get_bbox()
		{
			return m_bbox;
		}
		std::vector<float> get_points()
		{
			return m_points;
		}
		std::vector<float> get_segments()
		{
			return m_segments;
		}
		std::vector<float> get_triangulations()
		{
			return m_trigulations;
		}
		std::vector<Local_point> get_faces()
		{
			return m_faces;
		}
		std::vector<std::size_t> get_indexes()
		{
			return m_indexs;
		}
		std::vector<Local_point> get_all_points()
		{
			return m_all_points;
		}
	};
	//};
	template<class Nef_Polyhedron>
	ParseNef<Nef_Polyhedron>::ParseNef(Nef_Polyhedron& cv)
	{
		Volume_const_iterator c;
		Nef_Visitor SE(this);
		CGAL_forall_volumes(c, cv) {
			Shell_entry_const_iterator it;
			CGAL_forall_shells_of(it, c) {
				cv.visit_shell_objects(SFace_const_handle(it), SE);

			}
		}
	}
	template<class Nef_Polyhedron>
	ParseNef<Nef_Polyhedron>::ParseNef(std::vector<Nef_Polyhedron>& cv)
	{
		for (auto nef = cv.begin(); nef != cv.end(); ++nef)
		{
			Volume_const_iterator c;
			Nef_Visitor SE(this);
			CGAL_forall_volumes(c, *nef) {
				Shell_entry_const_iterator it;
				CGAL_forall_shells_of(it, c) {
					(*nef).visit_shell_objects(SFace_const_handle(it), SE);
				}
			}
		}
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::add_point(const Point& p1)
	{
		Local_point p = get_local_point(p1);
		m_all_points.push_back(p);
		m_bbox += p.bbox();
		m_points.push_back(static_cast<float>(p.x()));
		m_points.push_back(static_cast<float>(p.y()));
		m_points.push_back(static_cast<float>(p.z()));
		if (m_map.find(p) == m_map.end())
		{
			m_map[p] = index++;
		}
		assert(m_points.size() % 3 == 0);
	}
	template<class Nef_Polyhedron>
	inline void wxh::ParseNef<Nef_Polyhedron>::add_points(const Local_point & p1)
	{
		m_all_points.push_back(p1);
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::add_segment(const Point& p1, const Point& p2)
	{
		Local_point p11 = get_local_point(p1);
		Local_point p12 = get_local_point(p2);
		m_segments.push_back(static_cast<float>(p11.x()));
		m_segments.push_back(static_cast<float>(p11.y()));
		m_segments.push_back(static_cast<float>(p11.z()));
		m_segments.push_back(static_cast<float>(p12.x()));
		m_segments.push_back(static_cast<float>(p12.y()));
		m_segments.push_back(static_cast<float>(p12.z()));
		assert(m_segments.size() % 6 == 0);
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::add_point_in_triangulation(const Local_point& p1)
	{
		Local_point p = get_local_point(p1);
		m_trigulations.push_back(static_cast<float>(p.x()));
		m_trigulations.push_back(static_cast<float>(p.y()));
		m_trigulations.push_back(static_cast<float>(p.z()));
		if (m_map.find(p)==m_map.end())
		{
			m_map[p] = index++;
			m_indexs.push_back(m_map[p]);
		}
		else
		{
			m_indexs.push_back(m_map[p]);
		}
		assert(m_trigulations.size() % 3 == 0);
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::add_point_in_face(const Point& p1)
	{
		Local_point p = get_local_point(p1);
		m_faces.push_back(p);
		if (m_map.find(p) == m_map.end())
		{
			m_map[p] =index++;
		}
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::triangulate(Local_vector& normal, std::vector<Local_point>& vc)
	{
		if (vc.size()== 3)
		{
			triangular_face_end_internal(normal);
			return;
		}
		else if (is_facet_convex(vc, normal))
		{
			if (m_faces.size() == 4)
			{
				convex_quadrangular_face_end_internal(normal);
			}
			else
			{
				convex_face_end_internal(normal);
			}
		}
		else
		{
			nonconvex_face_end_internal(normal);
		}
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::triangular_face_end_internal(const Local_vector& normal)
	{
		for (int i = 0; i < 3; ++i)
		{
			add_point_in_triangulation(m_faces[i]);
		}
		//assert((index+1) % 3 == 0);
	}
	template<class Nef_Polyhedron>
	bool ParseNef<Nef_Polyhedron>::is_facet_convex(const std::vector<Local_point>& facet, const Local_vector& normal)
	{
		Local_kernel::Orientation orientation, local_orientation;
		std::size_t id = 0;
		do
		{
			const Local_point& S = facet[id];
			const Local_point& T = facet[(id + 1 == facet.size()) ? 0 : id + 1];
			Local_vector V1 = Local_vector((T - S).x(), (T - S).y(), (T - S).z());
			const Local_point& U = facet[(id + 2 >= facet.size()) ? id + 2 - facet.size() : id + 2];
			Local_vector V2 = Local_vector((U - T).x(), (U - T).y(), (U - T).z());

			orientation = Local_kernel::Orientation_3()(V1, V2, normal);
			// Is it possible that orientation==COPLANAR ? Maybe if V1 or V2 is very small ?
		} while (++id != facet.size() &&
			(orientation == CGAL::COPLANAR));

		//Here, all orientations were COPLANAR. Not sure this case is possible,
		// but we stop here.
		if (orientation == CGAL::COPLANAR)
		{
			return false;
		}

		// Now we compute convexness
		for (id = 0; id < facet.size(); ++id)
		{
			const Local_point& S = facet[id];
			const Local_point& T = facet[(id + 1 == facet.size()) ? 0 : id + 1];
			Local_vector V1 = Local_vector((T - S).x(), (T - S).y(), (T - S).z());

			const Local_point& U = facet[(id + 2 >= facet.size()) ? id + 2 - facet.size() : id + 2];
			Local_vector V2 = Local_vector((U - T).x(), (U - T).y(), (U - T).z());

			local_orientation = Local_kernel::Orientation_3()(V1, V2, normal);

			if (local_orientation != CGAL::ZERO)
			{
				if (local_orientation != orientation)
				{
					return false;
				}
			}
			else
			{
				if (CGAL::scalar_product(V1, V2) < 0)
				{
					return false;
				}  //TS and TU are opposite
			}
		}
		return true;
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::convex_quadrangular_face_end_internal(const Local_vector& normal)
	{
		add_point_in_triangulation(m_faces[0]);
		add_point_in_triangulation(m_faces[1]);
		add_point_in_triangulation(m_faces[2]);

		add_point_in_triangulation(m_faces[0]);
		add_point_in_triangulation(m_faces[2]);
		add_point_in_triangulation(m_faces[3]);
		//assert((index+1) % 3 == 0);
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::convex_face_end_internal(const Local_vector& normal)
	{
		for (std::size_t i = 1; i < m_faces.size() - 1; ++i)
		{
			Local_point& p0 = m_faces[0];
			Local_point& p1 = m_faces[i];
			Local_point& p2 = m_faces[i + 1];
			add_point_in_triangulation(p0);
			add_point_in_triangulation(p1);
			add_point_in_triangulation(p2);
			//assert((index+1) % 3 == 0);
		}
	}
	template<class Nef_Polyhedron>
	void ParseNef<Nef_Polyhedron>::nonconvex_face_end_internal(const Local_vector& normal)
	{
		typedef CGAL::Projection_traits_3<CGAL::Exact_predicates_inexact_constructions_kernel> P_traits;
		typedef CGAL::Triangulation_vertex_base_2< P_traits> Vb;
		//typedef CGAL::Triangulation_face_base_with_info_2<Face_info, P_traits>     Fb1;
		typedef CGAL::Constrained_triangulation_face_base_2<P_traits>         Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                        TDS;
		typedef CGAL::Exact_predicates_tag                                         Itag;
		typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>    CDT;
		std::vector<CDT::Constraint> vc_edges;
		P_traits cdt_traits(normal);
		CDT cdt(cdt_traits);
		//std:vector<typename CDT::Constraint> vc_edges;
		typedef CDT::Constraint Constrained;
		Local_point begin_p = m_faces[0];
		for (int i = 0; i < m_faces.size(); ++i)
		{
			if (i == m_faces.size() - 1)
			{
				vc_edges.push_back({ m_faces[i],begin_p });
			}
			else
			{
				vc_edges.push_back({ m_faces[i],m_faces[i + 1] });
			}
		}
		cdt.insert_constraints(vc_edges.begin(), vc_edges.end());
		for (auto&& it = cdt.finite_faces_begin(); it != cdt.finite_faces_end(); ++it)
		{
			for (int i = 0; i < 3; ++i)
			{
				add_point_in_triangulation(it->vertex(i)->point());
			}
		}
		//assert(m_indexs.size() % 3 == 0);
	}
	template<class Nef_Polyhedron>
	typename ParseNef<Nef_Polyhedron>::Local_vector ParseNef<Nef_Polyhedron>::get_vertex_normal(Vertex_const_handle vh)
	{
		Local_vector normal = CGAL::NULL_VECTOR;

		SHalfedge_const_iterator it = vh->shalfedges_begin();
		SHalfedge_const_handle end = it;
		do {
			Local_vector n = get_face_normal(it);
			normal = typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
			it = it->snext();
		} while (it != end);

		if (!typename Local_kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
		{
			normal = (typename Local_kernel::Construct_scaled_vector_3()(
				normal, 1.0 / CGAL::sqrt(normal.squared_length())));
		}

		return normal;
	}
	template<class Nef_Polyhedron>
	inline typename ParseNef<Nef_Polyhedron>::Local_vector wxh::ParseNef<Nef_Polyhedron>::get_face_normal(SHalfedge_const_handle she)
	{
		SHalfedge_around_facet_const_circulator he(she);
		Local_vector normal = CGAL::NULL_VECTOR;
		SHalfedge_around_facet_const_circulator end = he;
		unsigned int nb = 0;

		CGAL_For_all(he, end)
		{
			internal::newell_single_step_3(this->get_local_point
			(he->next()->source()->center_vertex()->point()),
				this->get_local_point(he->source()->center_vertex()->
					point()), normal);
			++nb;
		}

		assert(nb > 0);
		return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0 / nb));
	}
}
#endif