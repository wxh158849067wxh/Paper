#pragma once
#ifndef COMPUTE_H
#define COMPUTE_H
#include<vector>
#include<string>
#include<assert.h>
#include<unordered_map>
#include<unordered_set>
#include<utility>
#include<list>
#include<memory>
#include<set>
#include<unordered_set>

//#include <CGAL/Nef_polyhedron_3.h>
//#include<CGAL/Exact_predicates_exact_constructions_kernel.h>


#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h" 


	//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	//typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron_3;
namespace wxh {


	template<class Nef_polyhedron_3>
	class ComputeTO
	{
		typedef typename Nef_polyhedron_3::Vertex_const_handle Vertex_const_handle;
		typedef typename Nef_polyhedron_3::Halfedge_const_handle Halfedge_const_handle;
		typedef typename Nef_polyhedron_3::Halffacet_const_handle Halffacet_const_handle;
		typedef typename Nef_polyhedron_3::SHalfedge_const_handle SHalfedge_const_handle;
		typedef typename Nef_polyhedron_3::SHalfloop_const_handle SHalfloop_const_handle;
		typedef typename Nef_polyhedron_3::SFace_const_handle SFace_const_handle;
		typedef typename Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
		typedef typename Nef_polyhedron_3::Shell_entry_const_iterator Shell_entry_const_iterator;
		typedef typename Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
		typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
		typedef typename Nef_polyhedron_3::Halffacet_cycle_iterator Halffacet_cycle_iterator;
		typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator SHalfedge_around_facet_const_circulator;
		typedef typename Nef_polyhedron_3::Halffacet_const_iterator Halffacet_const_iterator;
		typedef typename Nef_polyhedron_3::SHalfedge_const_iterator SHalfedge_const_iterator;
		typedef typename Nef_polyhedron_3::Nef_polyhedron_S2 Nef_polyhedron_S2;
		typedef typename Nef_polyhedron_3::SHalfedge_handle SHalfedge_handle;
		typedef typename Nef_polyhedron_3::SHalfedge_around_sface_const_circulator SHalfedge_around_sface_const_circulator;
		typedef typename Nef_polyhedron_3::SFace_cycle_const_iterator SFace_cycle_const_iterator;
		typedef typename Nef_polyhedron_3::SVertex_const_iterator SVertex_const_iterator;
		typedef typename Nef_polyhedron_3::SHalfloop_const_iterator SHalfloop_const_iterator;
		typedef typename Nef_polyhedron_3::Shell_entry_iterator Shell_entry_iterator;
		typedef typename Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;
		typedef typename Nef_polyhedron_3::Point_3 Point;
		typedef typename Nef_polyhedron_3::Halfedge_const_iterator Halfedge_const_iterator;
	public:
		struct BasicNode {
			std::string type;//具体单元交类型
			std::vector<Point> data;//存储了这个基本交的数据
			short dimension;//表示维度 0表示点，1表示线，2表示面，3表示体
			int order = 0;//保留，可以表示交的序号
			Nef_polyhedron_3 bn_nef;//将每一个basicNode都包含一个nef;
			BasicNode() = default;
			BasicNode(const BasicNode&) = default;
			BasicNode& operator=(const BasicNode&) = default;
			BasicNode(std::string type1, std::vector<Point>& vc, short di)
				:type(std::move(type1)), data(std::move(vc)), dimension(di)
			{
			}
			bool operator==(const BasicNode& r) const
			{
				return this->type == r.type;
			}
			BasicNode(BasicNode&& node)
				:type(std::move(node.type)), data(std::move(node.data)), dimension(node.dimension), order(node.order),bn_nef(std::move(node.bn_nef))
			{}
			bool operator<(const BasicNode& r)
			{
				return this->order < r.order;
			}
			BasicNode& operator=(BasicNode&& node) = default;
		};
		struct hash_node
		{
			std::size_t operator()(const BasicNode& node) const
			{
				return std::hash<std::string>()(node.type);
			}
		};
		ComputeTO(Nef_polyhedron_3& nef1, Nef_polyhedron_3& nef2);//nef1 默认为体，nef2默认为线
		std::vector<BasicNode> compute();//返回结果
		~ComputeTO();
	private:
		Nef_polyhedron_3& m_nef1;
		Nef_polyhedron_3& m_nef2;
		std::shared_ptr< Nef_polyhedron_3>  m_union;//并的对象
		std::shared_ptr < Nef_polyhedron_3>  m_intersect_boundary;//交于外部的对象
		std::shared_ptr < Nef_polyhedron_3>  m_intersect_interior;//交于内部的对象
		std::shared_ptr < Nef_polyhedron_3>  m_intersect;//两者的交
		std::shared_ptr < Nef_polyhedron_3>  m_differ_lp;//线差体
		std::shared_ptr < Nef_polyhedron_3>  m_differ_pl;//体差线
		std::shared_ptr < Nef_polyhedron_3>  m_diff_LI;//线差线与体内相交的部分
		std::shared_ptr<Point> start_point;//nef2的起点
		std::shared_ptr<Point> end_point;//nef2的终点
		std::set<Point> m_set;//两个端点
		std::vector<BasicNode> result;//因为这是要返回的结果，此处没有用智能指针管理，避免类外不用智能指针导致内存出现问题，

		Vertex_const_iterator has_vertex(const Nef_polyhedron_3& nef, const  Vertex_const_iterator& vi);//重要函数，帮助判断vi是否在nef内，后续许多情况需要这个函数帮助判断
		Vertex_const_iterator has_vertex_p(const Nef_polyhedron_3& nef, const  Point& vi);//重要函数，帮助判断vi是否在nef内，后续许多情况需要这个函数帮助判断
		void find_boundary_points();//找到边体边界上孤立的基本交并且判断其类型
		void find_boundary_edges();//找到体边界上的线基本交并且判断其类型
		void find_interior_edges();//找到体内的线基本交并且判断其类型
		void find_endpoints();//找到线的端点；
		void get_nef(BasicNode& bn);
	};
	//nef1默认为体， nef默认为线
	template<class Nef_polyhedron_3>
	inline ComputeTO<Nef_polyhedron_3>::ComputeTO(Nef_polyhedron_3& nef1, Nef_polyhedron_3 &nef2)
		:m_nef1(nef1), m_nef2(nef2)
	{

		spdlog::get("basic_logger")->info("************************************************************************************");
		if ((m_nef1*m_nef2).is_empty())
		{
			spdlog::get("basic_logger")->info("type is NULL");
			return;
		}
		else
		{
			m_union = std::make_shared<Nef_polyhedron_3>(nef1 + nef2);
			m_intersect_boundary = std::make_shared<Nef_polyhedron_3>(nef1.boundary()*nef2);
			m_intersect_interior = std::make_shared<Nef_polyhedron_3>(nef1.interior()*nef2);
			m_intersect = std::make_shared<Nef_polyhedron_3>(nef1*nef2);
			m_differ_lp = std::make_shared<Nef_polyhedron_3>(nef2 - nef1);
			m_differ_pl = std::make_shared<Nef_polyhedron_3>(nef1 - nef2);
			m_diff_LI = std::make_shared<Nef_polyhedron_3>(nef2 - (*m_intersect_interior));
			find_endpoints();
			find_boundary_edges();
			find_boundary_points();
			find_interior_edges();
			//Nef_polyhedron_3 complete = (*m_intersect_boundary) + (*m_differ_lp) + (*m_intersect_interior);
			//assert(complete == nef2);
			std::map<Point, std::size_t> order_map;
			std::list<Point> order_list;
			int index = 1;
			for (auto it = nef2.vertices_begin(); it != nef2.vertices_end(); ++it)
			{
				order_list.push_back(it->point());
				order_map[it->point()] = index++;
			}
			for (auto it = result.begin(); it != result.end(); ++it)
			{
				if (it->type == "(2, 0, <null, out>)"||it->type== "(1, 1, <null, out>)"||it->type== "(2, 0, <on, out>)"||it->type== "(2, 0, <in, out>)")
				{
					Point& p1 = it->data.front();
					Point& p2 = it->data.back();
					if (order_map.count(p1) != 0 && order_map.count(p2) != 0)
					{
						continue;
					}
					Vertex_const_iterator vi1 = has_vertex_p(*m_differ_lp, p1);
					Vertex_const_iterator vi2 = has_vertex_p(*m_differ_lp, p2);
					//get_nef(*it);
					
					if (order_map.count(p2)==0&& order_map.count(p1)!=0)
					{
						assert(vi2->number_of_svertices() == 1);
						assert(vi1 == nullptr);
						Vertex_const_iterator tmp_vi = vi2->svertices_begin()->target();//得到被插入vertex的下一个
						Point p3 = tmp_vi->point();
						auto tmp_it = order_list.begin();
						for (; tmp_it != order_list.end(); ++tmp_it)
						{
							if ((*tmp_it) == p3)
							{
								break;
							}
						}
						assert(tmp_it != order_list.end());
						++tmp_it;//得到应该被插入的位置
						order_list.insert(tmp_it,p2);
					}
					else if (order_map.count(p1) == 0&& order_map.count(p2) != 0)
					{
						assert(vi1->number_of_svertices() == 1);
						assert(vi2 == nullptr);
						Vertex_const_iterator tmp_vi = vi1->svertices_begin()->target();//得到被插入vertex的下一个
						Point p3 = tmp_vi->point();
						auto tmp_it = order_list.begin();
						for (; tmp_it != order_list.end(); ++tmp_it)
						{
							if ((*tmp_it) == p3)
							{
								break;
							}
						}
						assert(tmp_it != order_list.end());
						++tmp_it;//得到应该被插入的位置
						order_list.insert(tmp_it,p1);
					}
					else
					{
						assert(false);
					}
				}
				else if (it->type =="(2, 0, <out, out>)" )
				{
					Point& p1 = it->data.front();
					Point& p2 = it->data.back();
					if (order_map.count(p1) != 0 && order_map.count(p2) != 0)
					{
						continue;
					}
					Vertex_const_iterator vi1 = has_vertex_p(*m_differ_lp, p1);
					Vertex_const_iterator vi2 = has_vertex_p(*m_differ_lp, p2);
					if (order_map.count(p1) == 0 && order_map.count(p2) == 0)
					{
						assert(vi1->number_of_svertices() == 1);
						assert(vi2->number_of_svertices() == 1);
						Vertex_const_iterator tmp_vi1 = vi1->svertices_begin()->target();//得到被插入vertex的下一个
						Point p31 = tmp_vi1->point();
						Vertex_const_iterator tmp_vi2 = vi2->svertices_begin()->target();//得到被插入vertex的下一个
						Point p32 = tmp_vi2->point();
						typename std::list<Point>::iterator tmp_it11 = std::list<Point>::iterator();
						typename std::list<Point>::iterator  tmp_it12 = std::list<Point>::iterator();
						bool mark1 = false;
						bool mark2 = false;
						auto tmp_it = order_list.begin();
						for (; tmp_it != order_list.end(); ++tmp_it)
						{
							if ((*tmp_it) == p31)
							{
								mark1 = true;
								tmp_it11 = tmp_it;
							}
							if ((*tmp_it) == p32)
							{
								mark2 = true;
								tmp_it12 = tmp_it;
							}
							if (mark1&&mark2)
							{
								break;
							}
						}
						assert(tmp_it11 != order_list.end());
						++tmp_it11;
						order_list.insert(tmp_it11,p1);
						++tmp_it12;
						order_list.insert(tmp_it12, p2);

					}
					else if (order_map.count(p1) == 0 && order_map.count(p2) != 0)
					{
						assert(vi1->number_of_svertices() == 1);
						assert(vi2->number_of_svertices() == 1);
						Vertex_const_iterator tmp_vi = vi1->svertices_begin()->target();//得到被插入vertex的下一个
						Point p3 = tmp_vi->point();
						auto tmp_it = order_list.begin();
						for (; tmp_it != order_list.end(); ++tmp_it)
						{
							if ((*tmp_it) == p3)
							{
								break;
							}
						}
						assert(tmp_it != order_list.end());
						++tmp_it;//得到应该被插入的位置
						order_list.insert(tmp_it, p1);
					}
					else if (order_map.count(p1) != 0 && order_map.count(p2) == 0)
					{
						assert(vi2->number_of_svertices() == 1);
						assert(vi1 == nullptr);
						Vertex_const_iterator tmp_vi = vi2->svertices_begin()->target();//得到被插入vertex的下一个
						Point p3 = tmp_vi->point();
						auto tmp_it = order_list.begin();
						for (; tmp_it != order_list.end(); ++tmp_it)
						{
							if ((*tmp_it) == p3)
							{
								break;
							}
						}
						assert(tmp_it != order_list.end());
						++tmp_it;//得到应该被插入的位置
						order_list.insert(tmp_it, p2);
					}
					else
					{
						assert(false);
					}
					
				}

			}
			order_map.clear();
			int index1 = 1;
			for (auto it = order_list.begin(); it != order_list.end(); ++it)
			{
				order_map[*it] = index1++;
			}
			int index2 = 1;
			for (auto vertex_h = order_list.begin(); vertex_h != order_list.end(); ++vertex_h)
			{
				//order_map[it->point()]=index++;
				typedef std::vector<Point>::iterator point_iterator;
				typedef std::pair<point_iterator, point_iterator> point_range;
				typedef std::list<point_range> polyline;
				//std::cout << *vertex_h<<" ";
				for (auto it = result.begin(); it != result.end(); ++it)
				{
					//std::cout << it->data.front() << std::endl;
					if (*vertex_h == (it->data).front())
					{
						if (it->dimension == 0)
						{
							if (it->bn_nef != Nef_polyhedron_3())
							{
								break;
							}
							it->order = index2++;
							get_nef(*it);
							break;
						}
						else if (it->dimension == 1)
						{

							if (it->bn_nef != Nef_polyhedron_3())
							{
								continue;
							}
							Point& start = it->data.front();
							Point& end = it->data.back();
							if (*vertex_h== start||*vertex_h==end)
							{
								it->order = index2++;
								get_nef(*it);
							}
						}
					}
				}

			}
			std::sort(result.begin(), result.end());
		}
	}
	template<class Nef_polyhedron_3>
	inline typename ComputeTO<Nef_polyhedron_3>::Vertex_const_iterator ComputeTO<Nef_polyhedron_3>::has_vertex(const Nef_polyhedron_3& nef, const  Vertex_const_iterator& vi)
	{
		for (Vertex_const_iterator it = nef.vertices_begin(); it != nef.vertices_end(); ++it)
		{
			if (it->point() == vi->point())
			{
				return it;
			}

		}
		return Vertex_const_iterator();
	}
	template<class Nef_polyhedron_3>
	inline  typename ComputeTO<Nef_polyhedron_3>::Vertex_const_iterator ComputeTO<Nef_polyhedron_3>::has_vertex_p(const Nef_polyhedron_3 & nef, const Point & vi)
	{
		for (Vertex_const_iterator it = nef.vertices_begin(); it != nef.vertices_end(); ++it)
		{
			if (it->point() == vi)
			{
				return it;
			}
		}
		return Vertex_const_iterator();
	}
	template<class Nef_polyhedron_3>
	inline void ComputeTO<Nef_polyhedron_3>::find_boundary_points()
	{
		for (auto it = m_intersect_boundary->vertices_begin(); it != m_intersect_boundary->vertices_end(); ++it)
		{
			if (it->number_of_svertices() == 0 && has_vertex(*m_intersect_interior, it) == nullptr)
			{//当前点为基本交
				Vertex_const_iterator p_out = has_vertex(*m_differ_lp, it);
				Vertex_const_iterator p_in = has_vertex(*m_intersect_interior, it);
				assert(p_in == nullptr);
				assert(p_out != nullptr);
				if (p_out->number_of_svertices() == 2)
				{
					assert(m_set.count((it)->point()) == 0);//该点不为端点
					std::vector<Point> vc(1, (it)->point());
					result.emplace_back("(1, 1, <out, out>)", vc, 0);
					spdlog::get("basic_logger")->info("体边界的点基本交是({0} {1} {2})", CGAL::to_double(it->point().x()), CGAL::to_double(it->point().y()), CGAL::to_double(it->point().z()));
					spdlog::get("basic_logger")->info("type is {}", result.back().type);
				}
				else if (p_out->number_of_svertices() == 1)
				{
					assert(m_set.count((it)->point()) == 1);//该点为端点
					std::vector<Point> vc(1, (it)->point());
					result.emplace_back("(1, 1, <out, null>)", vc, 0);
					spdlog::get("basic_logger")->info("体边界的点基本交是({0} {1} {2})", CGAL::to_double(it->point().x()), CGAL::to_double(it->point().y()), CGAL::to_double(it->point().z()));
					spdlog::get("basic_logger")->info("type is {}", result.back().type);
				}
				else
					assert(false);

			}

		}
	}
	template<class Nef_polyhedron_3>
	inline void ComputeTO<Nef_polyhedron_3>::find_boundary_edges()
	{
		std::unordered_set<Halfedge_const_iterator> u_set;
		for (auto it = m_intersect_boundary->halfedges_begin(); it != m_intersect_boundary->halfedges_end(); ++it)
		{
			if (it->source()->number_of_svertices() == 1)//只是加入端点
			{
				if (u_set.find(it) == u_set.end())//判断该端点是否已经处理过
				{
					Vertex_const_iterator vo = it->source();//线基本交的一端
					std::vector<Point> vc;
					vc.push_back(vo->point());
					auto tmp1 = it;
					Vertex_const_iterator vi = tmp1->target();//边界线基本交的另一端
					vc.push_back(vi->point());
					u_set.insert(tmp1);
					u_set.insert(tmp1->twin());
					while (vi->number_of_svertices() != 1)
					{
						auto reverse = tmp1->twin();
						tmp1 = (reverse == vi->svertices_begin() ? vi->svertices_last() : vi->svertices_begin());
						u_set.insert(tmp1);
						u_set.insert(tmp1->twin());
						vi = tmp1->target();
						vc.push_back(vi->point());
					}
					assert(vo->number_of_svertices() == 1);//以下的assert为保证边界线基本交正确
					assert(vi->number_of_svertices() == 1);
					assert(vo->point() == vc.front());
					assert(vi->point() == vc.back());
					assert(vc.size() >= 2);//now vc的两个端点就是线边界基本交的两个端点
					spdlog::get("basic_logger")->info("边界线基本交的两个端点分别是 1->({0} {1} {2}) ,  2->( {3} {4} {5} )", CGAL::to_double(vo->point().x()), CGAL::to_double(vo->point().y()), CGAL::to_double(vo->point().z()), CGAL::to_double(vi->point().x()), CGAL::to_double(vi->point().y()), CGAL::to_double(vi->point().z()));
					//vc.push_back(tmp1->point());
					//判断的思想主要是判断线基本交的两个端点是否属于线差体、是否属于线体的内部交、是否是端点
					auto res_out1 = has_vertex(*m_differ_lp, vo);//线和体的差
					auto res_out2 = has_vertex(*m_differ_lp, vi);//线个体的差
					auto res_in1 = has_vertex(*m_intersect_interior, vo);
					auto res_in2 = has_vertex(*m_intersect_interior, vi);
					bool is_out_1 = (res_out1 == nullptr ? false : true);
					bool is_out_2 = (res_out2 == nullptr ? false : true);
					bool is_in_1 = (res_in1 == nullptr ? false : true);
					bool is_in_2 = (res_in2 == nullptr ? false : true);
					bool is_endpoint1 = ((m_set.find(vo->point()) != m_set.end()));//是否为端点，true表示为端点
					bool is_endpoint2 = (m_set.find(vi->point()) != m_set.end());
					if (!is_out_1 && !is_out_2 && !is_in_1 && !is_in_2&&is_endpoint1&&is_endpoint2)
					{
						result.emplace_back("(2, 1, <null, null>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if ((!is_out_1&&is_out_2 && !is_in_1 && !is_in_2&&is_endpoint1 && !is_endpoint2) || (is_out_1 && !is_out_2 && !is_in_1 && !is_in_2 && !is_endpoint1&&is_endpoint2))
					{
						result.emplace_back("(2,1,<null,out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if ((!is_out_1 && !is_out_2 && !is_in_1&&is_in_2&&is_endpoint1 && !is_endpoint2) || (!is_out_1 && !is_out_2&&is_in_1 && !is_in_2 && !is_endpoint1&&is_endpoint2))
					{
						result.emplace_back("(2,1,<null,in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_out_1 && !is_out_2&&is_in_1&&is_in_2 && !is_endpoint1 && !is_endpoint2)
					{
						result.emplace_back("(2,1,<in,in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_out_1&&is_out_2 && !is_in_1 && !is_in_2 && !is_endpoint1 && !is_endpoint2)
					{
						result.emplace_back("(2,1,<out,out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_out_1&&is_out_2&&is_in_1 && !is_in_2 && !is_endpoint1 && !is_endpoint2 || is_out_1 && !is_out_2 && !is_in_1&&is_in_2 && !is_endpoint1 && !is_endpoint2)
					{
						result.emplace_back("(2,1,<in,out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else
						assert(false);
				}
			}
		}
	}
	template<class Nef_polyhedron_3>
	inline void ComputeTO<Nef_polyhedron_3>::find_interior_edges()
	{
		std::unordered_set< Halfedge_const_iterator> u_set;
		for (Halfedge_const_iterator it = m_intersect_interior->halfedges_begin(); it != m_intersect_interior->halfedges_end(); ++it)
		{
			if (it->source()->number_of_svertices() == 1 || it->source()->mark() == false)//只有当这条边的source的svertex个数为1或者，source的mark为false
			{
				if (u_set.find(it) == u_set.end())
				{
					std::vector<Point> vc;
					Vertex_const_iterator vl = it->source();
					vc.push_back(vl->point());
					u_set.insert(it);
					u_set.insert(it->twin());
					Halfedge_const_iterator h1 = it;
					Vertex_const_iterator vr = h1->target();
					vc.push_back(vr->point());
					while (vr->number_of_svertices() != 1 && vr->mark() == true)
					{
						Halfedge_const_iterator reverse = h1->twin();
						h1 = (reverse == vr->svertices_begin() ? vr->svertices_last() : vr->svertices_begin());
						u_set.insert(h1);
						u_set.insert(h1->twin());
						vr = h1->target();
						vc.push_back(vr->point());
					}
					//以下assert为确保基本交的正确性
					assert((vl->number_of_svertices() == 1 || vl->mark() == false));
					assert((vr->number_of_svertices() == 1 || vr->mark() == false));
					assert(vl->point() == vc.front() && vr->point() == vc.back());
					spdlog::get("basic_logger")->info("内部线基本交的两个端点分别是 1->({0} {1} {2}) ,  2->( {3} {4} {5} )", CGAL::to_double(vl->point().x()), CGAL::to_double(vl->point().y()), CGAL::to_double(vl->point().z()), CGAL::to_double(vr->point().x()), CGAL::to_double(vr->point().y()), CGAL::to_double(vr->point().z()));
					Vertex_const_iterator v_out_1 = has_vertex(*m_differ_lp, vl);
					Vertex_const_iterator v_out_2 = has_vertex(*m_differ_lp, vr);
					Vertex_const_iterator v_on_1 = has_vertex(*m_intersect_boundary, vl);
					Vertex_const_iterator v_on_2 = has_vertex(*m_intersect_boundary, vr);
					Vertex_const_iterator v_in_1 = has_vertex(*m_intersect_interior, vl);
					Vertex_const_iterator v_in_2 = has_vertex(*m_intersect_interior, vr);
					bool is_out1 = (v_out_1 == nullptr ? false : true);
					bool is_out2 = (v_out_2 == nullptr ? false : true);
					bool is_on1 = (v_on_1 == nullptr ? false : true);
					bool is_on2 = (v_on_2 == nullptr ? false : true);
					bool is_in1 = (v_in_1->number_of_svertices() == 1 ? false : true);//v_in_1或v_in_2肯定不为nullptr;如果number_of_svertices==1则说明该点没有邻接的体内单元交，反之==2则该点有邻接的体内单元交
					bool is_in2 = (v_in_2->number_of_svertices() == 1 ? false : true);//number_of_svertices肯定不为0
					assert(v_in_1->number_of_svertices() != 0);
					assert(v_in_2->number_of_svertices() != 0);

					bool is_endpoint1 = (m_set.find(vl->point()) != m_set.end());//是否为端点，是为true；
					bool is_endpoint2 = (m_set.find(vr->point()) != m_set.end());
					if (is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(0, 2, <null, null>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1&&is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2 || !is_endpoint1&&is_endpoint2&&is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <null, out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1&&is_in2 || !is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2&&is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <null, in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2 || !is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <null, on>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1&&is_in2 || !is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2&&is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <on, in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2&&is_out1&&is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <out, out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2 && !is_out1&&is_out2&&is_on1&&is_on2&&is_in1 && !is_in2 || !is_endpoint1 && !is_endpoint2&&is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1&&is_in2)
					{
						result.emplace_back("(2, 0, <in, out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2&&is_in1&&is_in2)
					{
						result.emplace_back("(2, 0, <in, in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2 && !is_out1&&is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2 || !is_endpoint1 && !is_endpoint2&&is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <on, out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (!is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2&&is_on1&&is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(2, 0, <on, on>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2 && !is_on1&&is_on2 && !is_in1 && !is_in2 || is_endpoint1 && is_endpoint2 && !is_out1 && !is_out2&&is_on1 && !is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(1, 1, <null, null>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1&&is_out2 && !is_on1&&is_on2 && !is_in1 && !is_in2 || !is_endpoint1&&is_endpoint2&&is_out1 && !is_out2&&is_on1 && !is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(1, 1, <null, out>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2 && !is_on1&&is_on2 && !is_in1&&is_in2 || !is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2&&is_on1 && !is_on2&&is_in1 && !is_in2)
					{
						result.emplace_back("(1, 1, <null, in>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1 && !is_endpoint2 && !is_out1 && !is_out2 && !is_on1&&is_on2 && !is_in1 && !is_in2 || !is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2&&is_on1 && !is_on2 && !is_in1 && !is_in2)
					{
						result.emplace_back("(1, 1, <null, on>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else if (is_endpoint1&&is_endpoint2 && !is_out1 && !is_out2 && !is_on1 && !is_on2 && !is_in1 && !is_in2 && (m_nef1 > m_nef2))
					{
						assert(m_nef1 > m_nef2);
						result.emplace_back("(0, 2, <null, null>)", vc, 1);
						spdlog::get("basic_logger")->info("type is {}", result.back().type);
					}
					else
						assert(false);




				}
			}
		}
	}

	template<class Nef_polyhedron_3>
	inline void ComputeTO<Nef_polyhedron_3>::find_endpoints()
	{
		for (auto it = m_nef2.vertices_begin(); it != m_nef2.vertices_end(); ++it)
		{
			assert(it->number_of_svertices() < 3);
			if (it->number_of_svertices() == 1)
			{
				m_set.insert(it->point());
			}
		}
	}

	template<class Nef_polyhedron_3>
	inline void ComputeTO<Nef_polyhedron_3>::get_nef( BasicNode& bn)
	{
		assert(bn.bn_nef == Nef_polyhedron_3());
		if (bn.dimension == 0)
		{
			assert(bn.data.size() == 1);
			bn.bn_nef = Nef_polyhedron_3(bn.data.front());
		}
		else if (bn.dimension == 1)
		{
			typedef typename std::vector<Point>::iterator point_iterator;
			typedef std::pair<point_iterator, point_iterator> point_range;
			typedef std::list<point_range> polyline;
			polyline line;
			line.push_back({ bn.data.begin(),bn.data.end() });
			bn.bn_nef = Nef_polyhedron_3(line.begin(), line.end(), Nef_polyhedron_3::Polylines_tag());
		}
		else
		{
			assert(false);
		}
	}

	template<class Nef_polyhedron_3>
	inline  std::vector<typename ComputeTO<Nef_polyhedron_3>::BasicNode> ComputeTO<Nef_polyhedron_3>::compute()
	{
		return result;
	}

	template<class Nef_polyhedron_3>
	inline ComputeTO<Nef_polyhedron_3>::~ComputeTO()
	{

	}


}
#endif // !COMPUTE_H

