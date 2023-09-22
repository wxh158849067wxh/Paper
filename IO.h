#pragma once
#ifndef IOFUNC_H
#define  IOFUNC_H
#include<string>
#include<fstream>
#include<sstream>
#include<utility>
#include<list>
#include<vector>
#include<assert.h>
#include<CGAL\Nef_nary_union_3.h>
#include<CGAL/OFF_to_nef_3.h>
#include <boost/filesystem.hpp>
namespace wxh {






	template<class Nef_polyhedron>
	std::vector<Nef_polyhedron> construct_multi_nef_from_dir(std::string dir_path)
	{
		boost::filesystem::path m_path(dir_path);
		if (!boost::filesystem::exists(m_path))
		{
			std::cerr << "does not exist!" << std::endl;
			assert(false);
		}
		std::list<Nef_polyhedron> res;
		/*if (!boost::filesystem::is_directory(m_path));
		{
			std::cerr << "not a path to dir" << std::endl;
			exit(-1);
		}*/
		for (auto&& x : boost::filesystem::directory_iterator(m_path))
		{
			Nef_polyhedron nef;
			boost::filesystem::path file_path = x.path();
			std::string file_name = file_path.filename().string();
			if (file_name.find(".off") != std::string::npos|| file_name.find(".OFF") != std::string::npos)
			{
				std::fstream off_file(file_path.string());
				if (!off_file.is_open())
				{
					std::cerr << "文件路径不存在，或出错" << std::endl;
				}
				CGAL::OFF_to_nef_3(off_file, nef);
				res.push_back(std::move(nef));
				off_file.close();
			}
			else if (file_name.find(".snc") != std::string::npos || file_name.find(".SNC") != std::string::npos)
			{
				std::fstream snc_file(file_path.string());
				if (!snc_file.is_open())
				{
					std::cerr << "文件路径不存在，或出错" << std::endl;
				}
				snc_file >> nef;
				res.push_back(std::move(nef));
				snc_file.close();
			}
			else if (file_name.find(".txt") != std::string::npos)
			{
				std::list<Nef_polyhedron> tmp_line= construct_multi_polyline_from_file<Nef_polyhedron>(file_name);
				res.assign(tmp_line.begin(), tmp_line.end());
			}
		}
		return res;
	}
	//从文件中构造一个可以有多个折线的nef polyhedron
	

	//从一个线段文件中构造多个Nef，一个折线为一个Nef，为接下来的计算做准备
	template<class Nef_polyhedron_3>
	std::list<Nef_polyhedron_3> construct_multi_polyline_from_file(std::string filename)
	{
		typedef typename Nef_polyhedron_3::Point_3 Point;
		typedef typename std::vector<Point>::iterator point_iterator;
		typedef std::pair<point_iterator, point_iterator> point_range;
		typedef std::list<point_range> polyline;
		std::list<Nef_polyhedron_3> res;
		std::fstream file;
		std::vector<Point> vc;
		file.open(filename);
		if (!file.is_open())
		{
			perror("file erro");
		}
		std::string line;
		int linenum = 0;
		polyline poly;
		while (getline(file, line))
		{
			if (line == ""&&linenum == 0)
			{
				continue;
			}
			++linenum;
			double x, y, z;
			if ((line == "" || line[0] == ' ') && linenum != 0)
			{
				assert(vc.size() >= 2);
				poly.push_back(point_range(vc.begin(), vc.end()));
				res.emplace_back(poly.begin(), poly.end(), typename Nef_polyhedron_3::Polylines_tag());
				vc.clear();
				poly.clear();
				linenum = 0;
				continue;
			}
			std::stringstream linedata(line);
			linedata >> x >> y >> z;
			vc.emplace_back(x, y, z);
		}
		file.close();
		if (vc.size() >= 2)
		{
			poly.push_back(point_range(vc.begin(), vc.end()));
			res.emplace_back(poly.begin(), poly.end(), typename Nef_polyhedron_3::Polylines_tag());
		}
		file.close();
		return res;
	}

	template<class Nef_polyhedron_3>
	Nef_polyhedron_3 construct_mutil_polyline_nef_from_file(std::string filename)
	{
		std::vector<Nef_polyhedron_3> vc = construct_multi_polyline_from_file<Nef_polyhedron_3 >(filename);
		CGAL::Nef_nary_union_3<Nef_polyhedron_3> nef_union;
		for (auto it = vc.begin(); it != vc.end(); ++it)
		{
			nef_union.add_polyhedron(std::move(*it));
		}
		return nef_union.get_union();
	}
	//从一个vector<BasicNode>的每个元素构造一个Nef polyhedron
	template <typename Nef_polyhedron_3>
	class  ComputeTO;
	template<class Nef_polyhedron_3>
	std::vector<Nef_polyhedron_3> construct_vec_of_nef_from_basic_node(std::vector < typename ComputeTO<Nef_polyhedron_3>::BasicNode>& b_vc)
	{
		typedef typename Nef_polyhedron_3::Point_3 Point;
		typedef typename std::vector<Point>::iterator point_iterator;
		typedef std::pair<point_iterator, point_iterator> point_range;
		typedef std::list<point_range> polyline;
		std::vector<Nef_polyhedron_3> res;
		for (auto it = b_vc.begin(); it != b_vc.end(); ++it)
		{
			polyline line;
			if (it->dimension == 0)
			{
				res.emplace_back((*it).data[0]);
			}
			else
			{
				line.push_back(point_range(it->data.begin(), it->data.end()));
				res.emplace_back(line.begin(), line.end(), typename Nef_polyhedron_3::Polylines_tag());
			}
		}
		return res;

	}

	template<class Nef_Polyhedron>
	Nef_Polyhedron construct_polyline_from_file(std::string s)
	{
		typedef typename Nef_Polyhedron::Point_3 Point;
		typedef typename std::vector<Point>::iterator point_iterator;
		typedef std::pair<point_iterator, point_iterator> point_range;
		typedef std::list<point_range> polyline;
		polyline poly;
		int number = 0;
		std::vector<std::shared_ptr<std::vector<Point>>>vec_of_vec;
		vec_of_vec.push_back(std::make_shared<std::vector<Point>>());
		Point pt;
		std::fstream file;
		file.open(s.c_str(), std::ios::in);
		if (!file.is_open())
		{
			std::cerr << "can not open the file!" << std::endl;
			file.close();
			Nef_Polyhedron N;
			return N;
		}
		std::string str1;
		std::string str2;
		std::string str3;
		while (!file.eof())
		{

			file >> str1;
			file >> str2;
			file >> str3;
			if (str1 == "-" || str2 == "-" || str3 == "-")
			{

				assert((*vec_of_vec[number]).size() >= 2);
				poly.push_back(point_range((*vec_of_vec[number]).begin(), (*vec_of_vec[number]).end()));
				++number;
				vec_of_vec.push_back(std::make_shared<std::vector<Point>>());

			}
			else if (str1 == "" || str2 == "" || str3 == "")
				break;
			else
			{

				double coordx = std::stod(str1);
				double coordy = std::stod(str2);
				double coordz = std::stod(str3);
				pt = Point(coordx, coordy, coordz);
				(*vec_of_vec[number]).push_back(pt);
			}

		}
		assert((*vec_of_vec[number]).size() >= 2);
		poly.push_back(point_range((*vec_of_vec[number]).begin(), (*vec_of_vec[number]).end()));
		Nef_Polyhedron N(poly.begin(), poly.end(), Nef_Polyhedron::Polylines_tag());
		file.close();
		return N;
	}

	//template<class Nef_polyhedron_3>
	//class ComputeTO;
	template<typename Nef_polyhedron_3, typename ComputeTO>
	Nef_polyhedron_3 construct_nef_from_basicnode(typename ComputeTO::BasicNode& node)
	{
		typedef typename Nef_polyhedron_3::Point_3 Point;
		typedef std::vector<Point>::iterator point_iterator;
		typedef std::pair<point_iterator, point_iterator> point_range;
		typedef std::list<point_range> polyline;
		if (node.dimension == 0)
		{
			return Nef_polyhedron_3(Point(node.data[0]));
		}
		else if (node.dimension == 1)
		{
			polyline poly;
			poly.push_back(point_range(node.data.begin(), node.data.end()));
			return Nef_polyhedron_3(poly.begin(), poly.end(), Nef_Polyhedron::Polylines_tag());
		}
	}

}

#endif



