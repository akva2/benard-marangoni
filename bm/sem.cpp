#include "sem.h"

using namespace mnl;
using namespace std;

void removeHydrostaticMode(basics::Matrix& res, const SEMInverseMassEvaluator& eval, int i)
{
  basics::Matrix* temp  = utilities::g_manager.clone(res);

  for( int j=0;j<res.cols();++j )
    for( int k=0;k<res.rows();++k )
      (*temp)[j][k] = Real(1);

  Real area = res.getDotter()(*temp,eval.getJacobian()[i]);

  basics::multPointwise(*temp,res,eval.getJacobian()[i]);
  Real mean = temp->sum()/area;

  res -= mean;

  utilities::g_manager.unlock(temp);
}

void removeHydrostaticMode(basics::matrixStack& res, const SEMInverseMassEvaluator& eval)
{
  basics::matrixStack* temp  = utilities::g_manager.clone(res);
  basics::matrixStack* temp2 = utilities::g_manager.clone(res);

  *temp = 1;

  eval.evaluateMass(*temp2,*temp);
  MPIDotter dotter(eval.m_rank,eval.m_size);
  Real area = dotter(*temp,*temp2);

  eval.evaluateMass(*temp2,res);
  Real mean = dotter(*temp,*temp2);
  mean /= area;

  res -= mean;

  utilities::g_manager.unlock(temp);
  utilities::g_manager.unlock(temp2);
}

string errorReport(basics::matrixStack& u, const basics::geometryStack& geometry,
    const basics::function2D& eu, const basics::Vector& grid, Real t, 
    const SEMLaplacianEvaluator* eval)
{
  basics::matrixStack* exact = utilities::g_manager.clone(u);
  basics::matrixStack* u2    = utilities::g_manager.clone(u);

  geometry.evaluate(*exact,grid,eu,t,eval->m_division[eval->m_rank].elements);
  u -= *exact;
  *u2 = u;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  if( eval->m_bc == SEMLaplacianEvaluator::PERIODIC )
    geometry.periodicDssum(*u2);
  else
    geometry.dssum(*u2,eval->m_rank,eval->m_size);
  MPIGeometryDotter<basics::geometryStack> dotter(geometry,	
      eval->m_division[eval->m_rank].elements,
      eval->m_rank,eval->m_size);
  Real L2U2 = dotter(u,*u2);
  *u2 = *exact;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  if( eval->m_bc == SEMLaplacianEvaluator::PERIODIC )
    geometry.periodicDssum(*u2);
  else
    geometry.dssum(*u2,eval->m_rank,eval->m_size);
  Real L2E = dotter(*exact,*u2);

  stringstream str;
  str << "error (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval ) {
    eval->evaluate(*u2,u);
    Real H1U2 = dotter(u,*u2);
    eval->evaluate(*u2,*exact);
    Real H1E = dotter(*exact,*u2);
    str << "error (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(exact);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(basics::matricesStack& u, const basics::geometryStack& geometry,
    const basics::function3D& eu, const basics::Vector& grid, Real t, 
    const extrudedLaplacianEvaluator* eval)
{
  basics::matricesStack* exact = utilities::g_manager.clone(u);
  basics::matricesStack* u2    = utilities::g_manager.clone(u);

  geometry.evaluate(*exact,grid,eu,t,eval->m_division[eval->m_rank].elements);
  u -= *exact;
  mask(u,geometry,eval->m_rank,eval->m_size,eval->m_bc);
  *u2 = u;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  dssum(*u2,geometry,eval->m_rank,eval->m_size,eval->m_bc);
  MPIGeometryDotter<basics::geometryStack> dotter(geometry,
      eval->m_division[eval->m_rank].elements,
      eval->m_rank,eval->m_size);
  Real L2U2 = dotter(u,*u2);
  *u2 = *exact;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  dssum(*u2,geometry,eval->m_rank,eval->m_size);
  Real L2E = dotter(*exact,*u2);

  stringstream str;
  str << "error (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval ) {
    eval->evaluate(*u2,u);
    Real H1U2 = dotter(u,*u2);
    eval->evaluate(*u2,*exact);
    Real H1E = dotter(*exact,*u2);
    str << "error (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(exact);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(basics::matricesStack& u, const basics::geometryStack3D& geometry,
    const basics::function3D& eu, const basics::Vector& grid, Real t, 
    const SEM3DLaplacianEvaluator* eval)
{
  basics::matricesStack* exact = utilities::g_manager.clone(u);
  basics::matricesStack* u2    = utilities::g_manager.clone(u);

  geometry.evaluate(*exact,eu,t,eval->m_division[eval->m_rank].elements);
  u -= *exact;
  *u2 = u;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  dssum(*u2,geometry,eval->m_rank,eval->m_size);
  MPIGeometryDotter<basics::geometryStack3D> dotter(geometry,
      eval->m_division[eval->m_rank].elements,
      eval->m_rank,eval->m_size);
  Real L2U2 = dotter(u,*u2);
  *u2 = *exact;
  geometry.mass(*u2,eval->m_division[eval->m_rank].elements);
  dssum(*u2,geometry,eval->m_rank,eval->m_size);
  Real L2E = dotter(*exact,*u2);

  stringstream str;
  str << "error (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval ) {
    eval->evaluate(*u2,u);
    Real H1U2 = dotter(u,*u2);
    eval->evaluate(*u2,*exact);
    Real H1E = dotter(*exact,*u2);
    str << "error (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(exact);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(mnl::basics::Field2<mnl::basics::matrixStack>& u,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function2D& ex, 
    const mnl::basics::function2D& ey,
    const mnl::basics::Vector& grid, Real t, 
    const SEMLaplacianEvaluator* eval2)
{
  basics::Field2<basics::matrixStack>* eu = utilities::g_manager.clone(u);
  basics::Field2<basics::matrixStack>* u2 = utilities::g_manager.clone(u);

  geometry.evaluate(*eu,grid,ex,ey,t,eval2->m_division[eval2->m_rank].elements);

  MPIGeometryDotter<basics::geometryStack> dotter(geometry,
      eval2->m_division[eval2->m_rank].elements,
      eval2->m_rank,eval2->m_size);
  u -= *eu;
  *u2 = u;
  geometry.mass(*u2,eval2->m_division[eval2->m_rank].elements);
  geometry.dssum(*u2,eval2->m_rank,eval2->m_size);
  Real L2U2  = dotter(u.X(),u2->X());
  L2U2 += dotter(u.Y(),u2->Y());

  *u2 = *eu;
  geometry.mass(*u2,eval2->m_division[eval2->m_rank].elements);
  geometry.dssum(*u2,eval2->m_rank,eval2->m_size);
  Real L2E  = dotter(eu->X(),u2->X());
  L2E += dotter(eu->Y(),u2->Y());

  stringstream str;
  str << "error velocity (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval2 ) {
    eval2->evaluateField(*u2,u);
    Real H1U2  = dotter(u.X(),u2->X());
    H1U2 += dotter(u.Y(),u2->Y());
    eval2->evaluateField(*u2,*eu);
    Real H1E   = dotter(eu->X(),u2->X());
    H1E  += dotter(eu->Y(),u2->Y());
    str << "error velocity (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(eu);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(mnl::basics::Field2<mnl::basics::matrixStack>& u,
    mnl::basics::matrixStack& p,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function2D& ex, 
    const mnl::basics::function2D& ey,
    const mnl::basics::function2D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const SEMInverseMassEvaluator& eval,
    const SEMLaplacianEvaluator* eval2)
{
  basics::matrixStack* pe = utilities::g_manager.clone(p);
  basics::matrixStack* p2 = utilities::g_manager.clone(p);

  geometry.evaluate(*pe,gridGL,ep,t, eval2->m_division[eval2->m_rank].elements);

  MPIDotter pdotter(eval2->m_rank,eval2->m_size);
  p -= *pe;
  removeHydrostaticMode(p,eval);
  eval.evaluateMass(*p2,p);
  Real L2P = pdotter(p,*p2);

  removeHydrostaticMode(*pe,eval);
  eval.evaluateMass(*p2,*pe);
  Real L2EP = pdotter(*pe,*p2);

  stringstream str;
  str << "error pressure (L^2): " << sqrt(L2P/L2EP) << endl;
  str << errorReport(u,geometry,ex,ey,grid,t,eval2);

  utilities::g_manager.unlock(pe);
  utilities::g_manager.unlock(p2);

  return( str.str() );
}

string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function3D& ex, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::Vector& grid, Real t, 
    const extrudedLaplacianEvaluator* eval2)
{
  basics::Field3<basics::matricesStack>* eu = utilities::g_manager.clone(u);
  basics::Field3<basics::matricesStack>* u2 = utilities::g_manager.clone(u);
  std::vector<int> desc = eval2->m_division[eval2->m_rank].elements;

  MPIGeometryDotter<basics::geometryStack> dotter(geometry,desc,eval2->m_rank,eval2->m_size);

  geometry.evaluate(*eu,grid,ex,ey,ez,t,desc);

  u -= *eu;
  //    u.save("hmms");
  *u2 = u;
  geometry.mass(*u2,desc);
  //    mask(*u2,geometry,0,1,eval2->m_bc);
  dssum(*u2,geometry,eval2->m_rank,eval2->m_size);
  Real L2U2  = dotter(u.X(),u2->X());
  L2U2 += dotter(u.Y(),u2->Y());
  L2U2 += dotter(u.Z(),u2->Z());

  *u2 = *eu;
  geometry.mass(*u2,desc);
  //    mask(*u2,geometry,0,1,eval2->m_bc);
  dssum(*u2,geometry,eval2->m_rank,eval2->m_size);
  Real L2E  = dotter(eu->X(),u2->X());
  L2E += dotter(eu->Y(),u2->Y());
  L2E += dotter(eu->Z(),u2->Z());

  stringstream str;
  str << "error velocity (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval2 ) {
    eval2->evaluateField(*u2,u);
    //        mask(*u2,geometry,0,1,eval2->m_bc);
    Real H1U2  = dotter(u.X(),u2->X());
    H1U2 += dotter(u.Y(),u2->Y());
    H1U2 += dotter(u.Z(),u2->Z());
    eval2->evaluateField(*u2,*eu);
    //        mask(*u2,geometry,0,1,eval2->m_bc);
    Real H1E   = dotter(eu->X(),u2->X());
    H1E  += dotter(eu->Y(),u2->Y());
    H1E  += dotter(eu->Z(),u2->Z());
    str << "error velocity (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(eu);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u,
    const mnl::basics::geometryStack3D& geometry,
    const mnl::basics::function3D& ex, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::Vector& grid, Real t, 
    const SEM3DLaplacianEvaluator* eval2)
{
  basics::Field3<basics::matricesStack>* eu = utilities::g_manager.clone(u);
  basics::Field3<basics::matricesStack>* u2 = utilities::g_manager.clone(u);
  std::vector<int> desc = eval2->m_division[eval2->m_rank].elements;

  MPIGeometryDotter<basics::geometryStack3D> dotter(geometry,desc,eval2->m_rank,eval2->m_size);

  geometry.evaluate(eu->X(),ex,t,desc);
  geometry.evaluate(eu->Y(),ey,t,desc);
  geometry.evaluate(eu->Z(),ez,t,desc);

  u -= *eu;
  //    u.save("hmms");
  *u2 = u;
  geometry.mass(*u2,desc);
  dssum(*u2,geometry,eval2->m_rank,eval2->m_size);
  Real L2U2  = dotter(u.X(),u2->X());
  L2U2 += dotter(u.Y(),u2->Y());
  L2U2 += dotter(u.Z(),u2->Z());

  *u2 = *eu;
  geometry.mass(*u2,desc);
  dssum(*u2,geometry,eval2->m_rank,eval2->m_size);
  Real L2E  = dotter(eu->X(),u2->X());
  L2E += dotter(eu->Y(),u2->Y());
  L2E += dotter(eu->Z(),u2->Z());

  stringstream str;
  str << "error velocity (L^2): " << sqrt(L2U2/L2E) << endl;

  if( eval2 ) {
    eval2->evaluateField(*u2,u);
    Real H1U2  = dotter(u.X(),u2->X());
    H1U2 += dotter(u.Y(),u2->Y());
    H1U2 += dotter(u.Z(),u2->Z());
    eval2->evaluateField(*u2,*eu);
    Real H1E   = dotter(eu->X(),u2->X());
    H1E  += dotter(eu->Y(),u2->Y());
    H1E  += dotter(eu->Z(),u2->Z());
    str << "error velocity (H^1): " << sqrt(H1U2/H1E) << endl;
  }

  utilities::g_manager.unlock(eu);
  utilities::g_manager.unlock(u2);

  return( str.str() );
}

string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u,
    mnl::basics::matricesStack& p,
    const mnl::basics::geometryStack& geometry,
    const mnl::basics::function3D& ex, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::function3D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const extrudedInverseMassEvaluator& eval,
    const extrudedLaplacianEvaluator* eval2)
{
  basics::matricesStack* pe = utilities::g_manager.clone(p);
  basics::matricesStack* p2 = utilities::g_manager.clone(p);
  std::vector<int> desc = eval2->m_division[eval2->m_rank].elements;

  MPIDotter pdotter(eval2->m_rank,eval2->m_size);

  geometry.evaluate(*pe,gridGL,ep,t,desc);
  p -= *pe;
  //    p.save("hmmp");
  removeHydrostaticMode(p,eval);
  eval.evaluateMass(*p2,p);
  Real L2P = pdotter(p,*p2);

  removeHydrostaticMode(*pe,eval);
  eval.evaluateMass(*p2,*pe);
  Real L2EP = pdotter(*pe,*p2);

  stringstream str;
  str << "error pressure (L^2): " << sqrt(L2P/L2EP) << endl;
  str << errorReport(u,geometry,ex,ey,ez,grid,t,eval2);

  utilities::g_manager.unlock(pe);
  utilities::g_manager.unlock(p2);

  return( str.str() );
}

string errorReport(mnl::basics::Field3<mnl::basics::matricesStack>& u,
    mnl::basics::matricesStack& p,
    const mnl::basics::geometryStack3D& geometry,
    const mnl::basics::function3D& ex, 
    const mnl::basics::function3D& ey,
    const mnl::basics::function3D& ez,
    const mnl::basics::function3D& ep,
    const mnl::basics::Vector& grid, 
    const mnl::basics::Vector& gridGL, Real t, 
    const SEM3DInverseMassEvaluator& eval,
    const SEM3DLaplacianEvaluator* eval2)
{
  basics::matricesStack* pe = utilities::g_manager.clone(p);
  basics::matricesStack* pt = utilities::g_manager.clone(u.X());
  basics::matricesStack* p2 = utilities::g_manager.clone(p);

  MPIDotter pdotter(eval2->m_rank,eval2->m_size);
  basics::Matrix GLL2G = utilities::GLL::interpolationMatrix(grid,gridGL);

  geometry.evaluate(*pt,ep,t);
  for( int i=0;i<pt->size();++i )
    geometry.interpolate((*pe)[i],(*pt)[i],GLL2G);

  p -= *pe;
  p -= p.sum()/p.length();
  eval.evaluateMass(*p2,p);
  Real L2P = pdotter(p,*p2);

  *pe -= pe->sum()/pe->length();
  eval.evaluateMass(*p2,*pe);
  Real L2EP = pdotter(*pe,*p2);

  stringstream str;
  str << "error pressure (L^2): " << sqrt(L2P/L2EP) << endl;
  str << errorReport(u,geometry,ex,ey,ez,grid,t,eval2);

  utilities::g_manager.unlock(pe);
  utilities::g_manager.unlock(p2);
  utilities::g_manager.unlock(pt);

  return( str.str() );
}

void G2GLL(basics::matricesStack& u, const basics::matricesStack& p,
    const basics::Vector& gridGLL, const basics::Vector& gridGL)
{
  int N = u[0].rows()-1;
  basics::Matrix G2GLL = utilities::GLL::interpolationMatrix(gridGL,gridGLL);
  basics::Matrix temp("yo",N+1,N-1);
  basics::Matrices temp2("yo",N+1,N+1,N-1);
  for( int i=0;i<p.size();++i ) {
    for( int l=0;l<N-1;++l ) {
      basics::multTranspose(temp,G2GLL,p[i][l],'N','N');
      basics::multTranspose(temp2[l],temp,G2GLL,'N','T');
    }
    basics::applyLocalGlobal(u[i],temp2,G2GLL,'N','T',0);
  }
}

Real doReduce(Real input, int rank, int size, int tag)
{
  Real result=input;
#ifdef HAS_MPI
  for( int i=0;i<size;++i ) {
    if( rank == i ) {
      for( int j=0;j<size;++j ) {
        if( rank == j )
          continue;
        MPI_Send(&input,1,MPI_DOUBLE,j,tag,MPI_COMM_WORLD);
      }
    } else {
      Real temp;
      MPI_Status res;
      MPI_Recv(&temp,1,MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&res);
      result += temp;
    }
  }
#endif

  return( result );
}

void vectorSum(basics::Vector& result, int rank, int size, int tag)
{
#ifdef HAS_MPI
  if( size == 1 )
    return;

  Real* temp  = new Real[result.length()];
  Real* temp2 = new Real[result.length()];
  memcpy(temp2,result.data(),result.length()*sizeof(Real));
  for( int i=0;i<size;++i ) {
    if( rank == i ) {
      for( int j=0;j<size;++j ) {
        if( rank == j )
          continue;
        MPI_Send(temp2,result.length(),MPI_DOUBLE,j,tag,MPI_COMM_WORLD);
      }
    } else {
      memset(temp,0,result.length()*sizeof(Real));
      MPI_Status res;
      MPI_Recv(temp,result.length(),MPI_DOUBLE,i,tag,MPI_COMM_WORLD,&res);
      result += temp;
    }
  }
  delete[] temp;
  delete[] temp2;
#endif
}

void sendStack(basics::matricesStack& result, const basics::matricesStack& input,
    const int* scount, const int* sdispl,
    const int* rcount, const int* rdispl,
    SEMLaplacianEvaluator::BC bc, vector<basics::matrixStack*>* stack,
    vector< vector<int> >* planes)
{
#ifdef HAS_MPI
  int skip=0;
  if(	bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM 	||
      bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    skip = 1;
  basics::Matrices* temp=NULL;
  if( planes ) {
    if( stack )
      temp = utilities::g_manager.clone(result[0]);
    else
      temp = utilities::g_manager.clone(input[0]);
  }

  for( int i=0;i<result.size();++i ) {
    if( stack ) {
      basics::Matrices* myresult = &result[i];
      if( planes )
        myresult = temp;
      MPI_Alltoallv(const_cast<Real*>(input[i][0].data()[0]),
          const_cast<int*>(scount),const_cast<int*>(sdispl),
          MPI_DOUBLE,(*myresult)[skip].data()[0],
          const_cast<int*>(rcount), const_cast<int*>(rdispl),
          MPI_DOUBLE,MPI_COMM_WORLD);
      if( planes ) {
        int l=skip;
        for( int j=0;j<planes->size();++j )
          for( int k=0;k<(*planes)[j].size();++k )
            result[i][(*planes)[j][k]+skip] = (*temp)[l++];
      }
    } else {
      basics::Matrices* myinput = const_cast<basics::Matrices*>(&input[i]);
      if( planes ) {
        myinput = temp;
        int l=skip;
        for( int j=0;j<planes->size();++j )
          for( int k=0;k<(*planes)[j].size();++k ) {
            (*temp)[l++] = input[i][(*planes)[j][k]+skip];
          }
      }
      MPI_Alltoallv((*myinput)[skip].data()[0],
          const_cast<int*>(scount), const_cast<int*>(sdispl),
          MPI_DOUBLE,result[i][0].data()[0],
          const_cast<int*>(rcount), const_cast<int*>(rdispl),
          MPI_DOUBLE,MPI_COMM_WORLD);
    }
  }

  if( stack ) {
    for( int l=0;l<stack->size();++l )
      delete (*stack)[l];
    stack->clear();
  }
  utilities::g_manager.unlock(temp);
#endif
}

void sendStack(basics::Field3<basics::matricesStack>& result,
    basics::matricesStack& buffer,
    const basics::matricesStack& input,
    const int* scount, const int* sdispl,
    const int* rcount, const int* rdispl,
    SEMLaplacianEvaluator::BC bc1,
    SEMLaplacianEvaluator::BC bc2,
    SEMLaplacianEvaluator::BC bc3,
    vector<basics::matrixStack*>* stack,
    vector< vector<int> >* planes)
{
  int len1=result[0][0].matrices();
  if(	bc1 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc1 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len1 -= 1;
  if( bc1 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len1 -= 2;
  int start1 = (len1-result[0][0].matrices())?1:0;
  int len2=result[0][0].matrices();
  if(	bc2 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc2 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len2 -= 1;
  if( bc2 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len2 -= 2;
  int start2 = (len2-result[1][0].matrices())?1:0;
  int len3=result[2][0].matrices();
  if(	bc3 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc3 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len3 -= 1;
  if( bc3 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len3 -= 2;
  int start3 = (len3-result[2][0].matrices())?1:0;

  sendStack(buffer,input,scount,sdispl,rcount,rdispl,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,stack,planes);

  for( int l=0;l<len1;++l )
    result[0].at(l+start1) = buffer.at(l);
  for( int l=0;l<len2;++l )
    result[1].at(l+start2) = buffer.at(l+len1);
  for( int l=0;l<len3;++l )
    result[2].at(l+start3) = buffer.at(l+len1+len2);
}

  vector<basics::matrixStack*>  
sendAndSetupStack(basics::matricesStack& result, const basics::matricesStack& input,
    SEMLaplacianEvaluator::BC bc, int rank, int size, 
    const int* scount, const int* sdispl, 
    const int* rcount, const int* rdispl,
    vector< vector<int> >* planes)
{
  int N = result[0].rows();

  sendStack(result,input,scount,sdispl,rcount,rdispl,bc,NULL,planes);
  vector<basics::matrixStack*> resultStack;
  for( int l=0;l<scount[rank]/(N*N);++l ) {
    vector<basics::Matrix*> foo;
    for( int i=0;i<size;++i ) {
      for( int k=0;k<result.size();++k )
        foo.push_back(&result[k][rdispl[i]/(N*N)+l]);
    }
    basics::matrixStack* temp = new basics::matrixStack("temp");
    temp->add(foo);
    resultStack.push_back(temp);
  }
  return( resultStack );
}

  vector<basics::matrixStack*>  
sendAndSetupStack(basics::matricesStack& result,
    basics::matricesStack& buffer,
    const basics::Field3<basics::matricesStack>& input,
    SEMLaplacianEvaluator::BC bc1, 
    SEMLaplacianEvaluator::BC bc2, 
    SEMLaplacianEvaluator::BC bc3, 
    int rank, int size, 
    const int* scount, const int* sdispl, 
    const int* rcount, const int* rdispl,
    vector< vector<int> >* planes)
{
  int len1=input[0][0].matrices();
  if(	bc1 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc1 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len1 -= 1;
  if( bc1 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len1 -= 2;
  int start1 = (input[0][0].matrices()-len1)?1:0;
  int len2=input[0][0].matrices();
  if(	bc2 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc2 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len2 -= 1;
  if( bc2 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len2 -= 2;
  int start2 = (input[1][0].matrices()-len2)?1:0;
  int len3=input[2][0].matrices();
  if(	bc3 == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM 		||
      bc3 == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM)
    len3 -= 1;
  if( bc3 == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    len3 -= 2;
  int start3 = (input[2][0].matrices()-len3)?1:0;

  for( int l=0;l<len1;++l )
    buffer.at(l) = input[0].at(l+start1);
  for( int l=0;l<len2;++l )
    buffer.at(l+len1) = input[1].at(l+start2);
  for( int l=0;l<len3;++l )
    buffer.at(l+len1+len2) = input[2].at(l+start3);

  return sendAndSetupStack(result,buffer,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,
      rank,size,scount,sdispl,rcount,rdispl,planes);
}

void massReference(basics::Matrix& res, const basics::Vector& weight)
{
  for( int j=0;j<res.cols();++j )
    res[j] *= weight[j];
  for( int k=0;k<res.rows();++k )
    res.scaleRow(k,weight[k]);
}

void massReference(basics::Matrices& res, const basics::Vector& weight)
{
  for( int l=0;l<res.matrices();++l ) {
    massReference(res[l],weight);
    res[l] *= weight[l];
  }
}

void dssum(mnl::basics::matricesStack& result, const mnl::basics::geometryStack& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc)
{
  if( size==1 ) {
    if( bc == SEMLaplacianEvaluator::PERIODIC ||
        bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET ||
        bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM )
      geometry.periodicDssum(result);
    else
      geometry.dssum(result);
    return;
  }

  std::vector<int*> scount, sdispl, rcount, rdispl;
  mnl::utilities::iterationStatistics stat(result[0].matrices(),size,
      scount,rcount,sdispl,rdispl);
  stat.getDisplacements(scount,sdispl,rcount,rdispl,rank,size,result[0][0].length());

  mnl::basics::matricesStack* temp = 
    mnl::utilities::g_manager.aquireMatricesStack("temp",result[0].rows(),result[0].cols(),
        (rdispl[0][size-1]+rcount[0][size-1])/result[0][0].length(),result.size());
  std::vector<mnl::basics::matrixStack*> stack = 
    sendAndSetupStack(*temp,result,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,
        rank,size,scount[0],sdispl[0],rcount[0],rdispl[0]);

  for( int l=0;l<stack.size();++l ) 
    geometry.dssum(*stack[l]);
  sendStack(result,*temp,scount[1],sdispl[1],rcount[1],rdispl[1],
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,&stack);
  stat.cleanDisplacements(scount,sdispl,rcount,rdispl);

  mnl::utilities::g_manager.unlock(temp);
}

void dssumCoarse(mnl::basics::matricesStack& result, int N,
    const mnl::basics::geometryStack& geometry, 
    int rank, int size, bool neumann)
{
  vector<mnl::basics::coarseGrid> group = geometry.getCoarseGroups(1);
  vector<mnl::basics::coarseGrid> group2 = geometry.getCoarseGroups(size);
  if( size==1 ) {
    geometry.dssumCoarse(result,group[0]);
    return;
  }

  if( rank == 0 ) {
    basics::matricesStack* stack = geometry.getCoarseBuffer(N,N,result[0].matrices(),group[0],
        "coarsetemp",false,neumann);
    for( unsigned int l=0;l<group2.size();++l ) {
      for( unsigned int k=0;k<group2[l].size();++k ) {
        if( l == 0 ) {
          (*stack)[group2[l][k].size3] = result[k];
        } else {
#ifdef HAS_MPI
          MPI_Recv((*stack)[group2[l][k].size3].data(),
              (*stack)[group2[l][k].size3].length(),
              MPI_DOUBLE,l,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
#endif
        }
      }
    }
    geometry.dssumCoarse(*stack,group[0],0,1,0);
    for( unsigned int l=0;l<group2.size();++l ) {
      for( unsigned int k=0;k<group2[l].size();++k ) {
        if( l == 0 ) {
          result[k] = (*stack)[group2[l][k].size3];
        }
        else {
#ifdef HAS_MPI
          MPI_Send((*stack)[group2[l][k].size3].data(),
              (*stack)[group2[l][k].size3].length(),
              MPI_DOUBLE,l,0,MPI_COMM_WORLD);
#endif
        }
      }
    }
    delete stack;
  } else {
    int sends=group2[rank].size();
#ifdef HAS_MPI
    for( int i=0;i<sends;++i )
      MPI_Send(result[i].data(),result[i].length(),MPI_DOUBLE,0,0,
          MPI_COMM_WORLD);
    for( int i=0;i<sends;++i ) {
      MPI_Recv(result[i].data(),result[i].length(),MPI_DOUBLE,0,0,
          MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif
  }
}

void dssum(mnl::basics::matricesStack& result, const mnl::basics::geometryStack3D& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc)
{
  if( size==1 ) {
    geometry.dssum(result);
    return;
  }

  mnl::basics::matricesStack* temp2 = geometry.localToGlobalZ(result,rank,size,false,true);

  std::vector<int*> scount, sdispl, rcount, rdispl;
  mnl::utilities::iterationStatistics stat((*temp2)[0].matrices(),size,
      scount,rcount,sdispl,rdispl);
  stat.getDisplacements(scount,sdispl,rcount,rdispl,rank,size,(*temp2)[0][0].length());


  mnl::basics::matricesStack* temp = 
    mnl::utilities::g_manager.aquireMatricesStack("temp",result[0].rows(),result[0].cols(),
        (rdispl[0][size-1]+rcount[0][size-1])/(*temp2)[0][0].length(),temp2->size());
  std::vector<mnl::basics::matrixStack*> stack = 
    sendAndSetupStack(*temp,*temp2,SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,
        rank,size,scount[0],sdispl[0],rcount[0],rdispl[0]);

  for( int l=0;l<stack.size();++l ) 
    geometry.dssum(*stack[l]);
  sendStack(*temp2,*temp,scount[1],sdispl[1],rcount[1],rdispl[1],
      SEMLaplacianEvaluator::HOMOGENOUS_NEUMANN,&stack);
  stat.cleanDisplacements(scount,sdispl,rcount,rdispl);

  mnl::utilities::g_manager.unlock(temp);
  geometry.globalToLocalZ(result,temp2,rank,size);
}

void mask(basics::matricesStack& result, const basics::geometryStack& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc)
{
  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    geometry.mask(result,rank,size);
  else if( bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM ) {
    result.at(0).clear();
    for( int l=1;l<result[0].matrices();++l )
      geometry.mask(result.at(l),rank,size);
  }
  else if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM )
    result.at(0).clear();
  else if( bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET ) {
    result.at(0).clear();
    result.at(result[0].matrices()-1).clear();
  } else if( bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM )
    result.at(0).clear();
}

void mask(basics::matricesStack& result, const basics::geometryStack3D& geometry,
    int rank, int size, SEMLaplacianEvaluator::BC bc)
{
  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET && geometry.m_level == 1) {
    geometry.mask(result,rank,size);
    return;
  }

  basics::coarseGrid desc = geometry.getDivisionInfo(size);

  vector<basics::matricesStack*> op = geometry.setupFakeStacks(result,geometry.m_level,desc[rank].elements.size()/geometry.m_level);

  int N = (*op[0])[0].matrices();

  if( bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM ||
      bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET ) {
    op[0]->at(0).clear();
    for( int i=0;i<op.size();++i )
      for( int l=0;l<(*op[i])[0].matrices();++l )
        geometry.mask(op[i]->at(l),rank,size);
  }
  else if( bc == SEMLaplacianEvaluator::NEUMANN_DIRICHLET_BOTTOM )
    op[0]->at(0).clear();
  if( bc == SEMLaplacianEvaluator::HOMOGENOUS_DIRICHLET )
    op[op.size()-1]->at(N-1).clear();

  geometry.killFakeStacks(op);
}

void mask(mnl::basics::Field3<mnl::basics::matricesStack>& result, 
    const basics::geometryStack& geometry, int rank, int size,
    SEMLaplacianEvaluator::BC bc)
{
  mask(result.X(),geometry,rank,size,bc);
  mask(result.Y(),geometry,rank,size,bc);
  mask(result.Z(),geometry,rank,size,bc);
  if( bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM ||
      bc == SEMLaplacianEvaluator::PERIODIC_DIRICHLET_BOTTOM )
    result.Z().at(result.Z()[0].matrices()-1).clear();
}

void mask(mnl::basics::Field3<mnl::basics::matricesStack>& result, 
    const basics::geometryStack3D& geometry, int rank, int size,
    SEMLaplacianEvaluator::BC bc)
{
  mask(result.X(),geometry,rank,size,bc);
  mask(result.Y(),geometry,rank,size,bc);
  mask(result.Z(),geometry,rank,size,bc);
  if( bc == SEMLaplacianEvaluator::DIRICHLET_DIRICHLET_BOTTOM) {
    basics::coarseGrid desc = geometry.getDivisionInfo(size);
    vector<basics::matricesStack*> op = geometry.setupFakeStacks(result.Z(),geometry.m_level,desc[rank].elements.size()/geometry.m_level);
    op[op.size()-1]->at((*op[0])[0].matrices()-1).clear();
    geometry.killFakeStacks(op);
  }
}
