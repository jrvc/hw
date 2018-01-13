/*
 * ex6.cpp
 *
 *  Created on: Jun 8, 2016
 *      Author: kahle
 */



//#include <bits/shared_ptr_base.h>
//#include <dolfin/adaptivity/ErrorControl.h>
#include <dolfin/common/Array.h>
//#include <dolfin/common/MPI.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/DirichletBC.h>
//#include <dolfin/fem/Equation.h>
#include <dolfin/fem/LinearVariationalSolver.h>
//#include <dolfin/fem/solve.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/generation/UnitSquareMesh.h>
#include <dolfin/io/File.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
//#include <dolfin/la/GenericVector.h>
#include <dolfin/log/log.h>
//#include <dolfin/mesh/Mesh.h>
//#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/la/solve.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/parameter/Parameter.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/refinement/PlazaRefinementND.h>
#include <dolfin/refinement/RegularCutRefinement.h>
#include <dolfin/refinement/UniformMeshRefinement.h>
#include <dolfin/la/PETScKrylovSolver.h>
#include <sys/stat.h>
#include <cassert>
#include <cmath>
#include <memory>
//#include <string>


#include "ex6_sys.h"


using std::endl;
using std::cout;




/*
 * this class describes the volume force f on the right hand side of the equation
 */
class VolumeForce : public dolfin::Expression
{

public:
	VolumeForce(unsigned geometric_dimension)
:Expression(),gdim(geometric_dimension){}

	void eval(dolfin::Array<double>& val, const dolfin::Array<double>& x)const
	{
		assert(val.size() == 1);
		assert(x.size() == gdim);
		assert(gdim==2);

		double x1 = x[0];
		double x2 = x[1];

		//using this rhs you will obtain an unexpected behaviour
		//		val[0] = 5*M_PI*M_PI*std::sin(M_PI*x1)*std::sin(2*M_PI*x2);

		val[0] = 1.0;


	}

private:
	const unsigned gdim;


};


/*
 * this class describes the Dirichlet boundary data u==g on \partial \Omega
 */
class BoundaryData : public dolfin::Expression
{

public:
	BoundaryData(unsigned geometric_dimension)
:Expression(),gdim(geometric_dimension){}

	void eval(dolfin::Array<double>& val, const dolfin::Array<double>& x)const
	{
		assert(val.size() == 1);
		assert(x.size() == gdim);
		assert(gdim==2);

		double x1 = x[0];
		double x2 = x[1];

		val[0] = 0.0;

	}

private:
	const unsigned gdim;


};



/*
 * this class describes the points on the boundary of \Omega, i.e. where to apply Dirichlet bonudary data
 */
class BoundaryDomain : public dolfin::SubDomain
{

public:

	bool inside(const dolfin::Array<double>& x, bool on_boundary)const
	{
		//every point x on the boundary is inside the Dirichlet part of the boundary
		return on_boundary;
	}

};




/*
 * create A and b, i.e. the finite dimensional representation of the PDE
 */
void assembleProblem(std::shared_ptr<dolfin::GenericMatrix>& A,
		std::shared_ptr<dolfin::GenericVector>& b,
		const std::shared_ptr<dolfin::FunctionSpace>& V)
{

	static VolumeForce force(V->mesh()->geometry().dim());
	static BoundaryData g(V->mesh()->geometry().dim());
	static BoundaryDomain bnd;

	A.reset(new dolfin::Matrix());
	b.reset(new dolfin::Vector());



	ex6_sys::BilinearForm a(V,V);
	dolfin::assemble(*A,a);




	ex6_sys::LinearForm f(V,force);
	dolfin::assemble(*b,f);

	dolfin::DirichletBC bc(V,dolfin::reference_to_no_delete_pointer(g),dolfin::reference_to_no_delete_pointer(bnd));

	bc.apply(*A,*b);


}


/*
 * overwrite mesh with a refined version
 */
void refineMesh(dolfin::Mesh& mesh)
{
	dolfin::Mesh mesh_old(mesh);
	dolfin::PlazaRefinementND::refine(mesh,mesh_old,true,false);

}




struct solver_cg_information
{
	unsigned number_of_iterations=0;
};


/*
 * EXERCISE: implement this function
 */
solver_cg_information solveLinearSystem_conjugatedGradients(
		const dolfin::GenericMatrix& A,dolfin::GenericVector& x,const dolfin::GenericVector& b,
		double tol)
{

	//guard to prevent from never stopping loop
	const unsigned maxit = x.size();



	unsigned m=0;//iteration count



	//initialize the iteration vectors,
	dolfin::Vector p_m; A.init_vector(p_m,0);
	dolfin::Vector r_m;  A.init_vector(r_m,0);
	dolfin::Vector v_m; A.init_vector(v_m,0);


	//'rename' x to x_m for the iteration
	dolfin::GenericVector& x_m = x;

	double lambda_m=1e16;


	A.mult(x_m,r_m);
	r_m.axpy(-1.0,b);
	r_m*=-1.0;

	p_m = r_m;

	double alpha_m = 1e16;
	double alpha_mp1 = r_m.inner(r_m);
	const double residuum_initial = std::sqrt(alpha_mp1);

	solver_cg_information res;


	dolfin::info("%.5s %16s","iter","residuum");

	for (m = 0;m<maxit;++m)
	{
		alpha_m = alpha_mp1;

		dolfin::info("%5d %16.14f",m,std::sqrt(alpha_m));

		if (std::sqrt(alpha_m) < tol*residuum_initial)
		{
			break;
		}

		//EXERCISE: fill the body of the algorithm.
		//NOTE: you will commonly need the y.axpy(a,x) operation, which means y:= a*x+y (for x,y Vectors, a Scalar)
		A.mult(p_m,v_m);
		
		// check this scalar product!!
		double scalProd = 0.;
		for (int i = 0; i < v_m.size(); ++i)
		{
		  scalProd = scalProd + v_m[i] * p_m[i]; 
		}
		lambda_m = alpha_m / scalProd;
		
		x_m.axpy(lambda_m, p_m);
		r_m.axpy(-lambda_m, v_m);
		alpha_mp1 = r_m.inner(r_m);
		p_m.axpy(alpha_mp1/alpha_m - 1,p_m);
		p_m.axpy(1,r_m);
		
	};



	res.number_of_iterations = m;



	return res;
}




solver_cg_information solveLinearSystem_dolfin(const dolfin::GenericMatrix& A,dolfin::GenericVector& x,const dolfin::GenericVector& b,
		double tol)
{

	/*
	 * some possible preconditioners:
	 * none,
	 * amg,
	 * sor
	 */

	//EXERCISE: choose the right preconditioner
	const std::string precond = "amg";

	dolfin::KrylovSolver sol("cg",precond);


	auto p = dolfin::KrylovSolver::default_parameters();

	//	dolfin::info(p,true);

	p["relative_tolerance"] = tol;
	p["absolute_tolerance"] = 0.0;
	p["maximum_iterations"] = (int)x.size();
	p["report"] = true;
	p["monitor_convergence"] = true;

	sol.update_parameters(p);




	solver_cg_information res;


	res.number_of_iterations = sol.solve(A,x,b);



	return res;

}





////////////////////////////////////////////////////////////////////////


/* measure the difference || u_h1-u_h2||_{L^2(\Omega)} ,  u_h2 is on the finer grid */
double difference_L2(const dolfin::Function& u1, const dolfin::Function& u2)
{


	//use finer solution for the form
	ex6_sys::Form_diff_L2 a(u1.function_space()->mesh());
	a.set_coefficient("u1",dolfin::reference_to_no_delete_pointer(u1));
	a.set_coefficient("u2",dolfin::reference_to_no_delete_pointer(u2));


	const double res_squared = dolfin::assemble(a);


	return  (res_squared > 0) ? std::sqrt(res_squared) : 0.0 ;

}

/* measure the difference || \nabla( u_h1-u_h2)||_{L^2(\Omega)} , u_h2 is on the finer grid */
double difference_H10(const dolfin::Function& u1, const dolfin::Function& u2)
{
	//use finer solution for the form
	ex6_sys::Form_diff_H10 a(u1.function_space()->mesh());
	a.set_coefficient("u1",dolfin::reference_to_no_delete_pointer(u1));
	a.set_coefficient("u2",dolfin::reference_to_no_delete_pointer(u2));


	const double res_squared = dolfin::assemble(a);
	return  (res_squared > 0) ? std::sqrt(res_squared) : 0.0 ;
}












/*
 * the main is the entry to the program
 */
int main(int argc, char* argv[])
{

	mkdir("output",0700);

	//maximum number of mesh refinements//max=8
	const unsigned max_lv = 7;


	//we solve the linear equation until the initial residuum is reduced by a factor of tol
	const double tol = 1e-12;


	//whether or not use investigate the self written cg code
	//EXERCISE: choose the corresponding value for the exercises
	const bool use_self_written_conjugated_gradients = false;

	//mesh and function space
	std::shared_ptr<dolfin::Mesh> mesh(new dolfin::UnitSquareMesh(4,4,"crossed"));
	std::shared_ptr<dolfin::FunctionSpace> V;


	//the solution using a direct solver (used as reference)
	std::shared_ptr<dolfin::Function> x_lu;
	//the solution using the iterative solver
	std::shared_ptr<dolfin::Function> x_cg;


	//the finite dimensional representation of the PDE
	std::shared_ptr<dolfin::GenericMatrix> A;
	std::shared_ptr<dolfin::GenericVector> b;

	dolfin::File solution_cg("output/uh_cg.pvd");




	dolfin::Array<unsigned> nIt(max_lv);
	dolfin::Array<unsigned> nDof(max_lv);

	for (unsigned lv = 0; lv < max_lv;++lv)
	{
		if (lv)
			refineMesh(*mesh);

		V.reset(new ex6_sys::FunctionSpace(mesh));

		dolfin::info("\nStart solving on a lv %u, with %u degrees of freedom",lv,V->dim());

		assembleProblem(A,b,V);
//now Ax=b is the finite dimensional representation

		x_cg.reset(new dolfin::Function(V));
		x_cg->vector()->zero();

		solver_cg_information info;
		if (use_self_written_conjugated_gradients)
			info = solveLinearSystem_conjugatedGradients(*A,*x_cg->vector(),*b,tol);
		else
			info = solveLinearSystem_dolfin(*A,*x_cg->vector(),*b,tol);

		nIt[lv] = info.number_of_iterations;
		nDof[lv] = V->dim();

  solution_cg<<*x_cg;

		if (V->dim()<10000)
		{
			x_lu.reset(new dolfin::Function(V));

			dolfin::solve(*A,*x_lu->vector(),*b,"lu","none");

			dolfin::info("||x_lu - x_cg||_{L2} = %g, ||x_lu - x_cg||_{H10} = %g \n\n",difference_L2(*x_lu,*x_cg),difference_H10(*x_lu,*x_cg));

			

		}

	}//end for all refinement level



	///

	dolfin::Array<double> rate(max_lv);
	rate[0] = 1;
	for (unsigned lv=1;lv<max_lv;++lv)
	{
		//EXERCISE: calculate the rate 'alpha' (NOTE: when you divide two numbers of type unsigned the result is an unsigned)
	  
	  
	}


	dolfin::info("%8s %8s %6s","ndof","nIt","alpha");
	for (unsigned lv = 0;lv<max_lv;++lv)
	{
		dolfin::info("%8u %8u %6.4f",nDof[lv],nIt[lv],rate[lv]);
	}




}


