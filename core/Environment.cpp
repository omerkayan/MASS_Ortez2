#include "Environment.h"
#include "DARTHelper.h"
#include "Character.h"
#include "BVH.h"
#include "Muscle.h"
#include "dart/collision/bullet/bullet.hpp"
using namespace dart;
using namespace dart::simulation;
using namespace dart::dynamics;
using namespace dart::math;
using namespace MASS;
/*
*library
std::ofstream tibia_tal_act;

*initialize
femur_tibia_aci.open ("femur_tibia_aci.txt", std::ios::app);

*reset bool rsi
tibia_femur_aci.close();
double x_tmp2;

*step
x_tmp2 = *(mCharacter->GetSkeleton()->getBodyNode("FemurL")->getChildJoint(0)->getPositions().data());
femur_tibia_aci<<x_tmp2<<std::endl;  
*/

std::ofstream reward_file;
std::ofstream activationlevel_file;
std::ofstream tibia_talus_aci;
std::ofstream femur_tibia_aci;
std::ofstream tibia_yer_aci;
std::ofstream head_neck_aci;
std::ofstream handR_acisal_hiz;
std::ofstream handL_acisal_hiz;

std::ofstream FemurL_x_desired_torque;
std::ofstream FemurL_y_desired_torque;
std::ofstream FemurL_z_desired_torque;
std::ofstream TibiaL_desired_torque;

Environment::
Environment()
	:mControlHz(30),mSimulationHz(900),mWorld(std::make_shared<World>()),mUseMuscle(true),w_q(0.65),w_v(0.1),w_ee(0.15),w_com(0.1)
{

}

std::shared_ptr<ArrowShape> mArrow; //ok kütüphanesi

void
Environment::
Initialize(const std::string& meta_file,bool load_obj)
{
	std::ifstream ifs(meta_file);
	if(!(ifs.is_open()))
	{
		std::cout<<"Can't read file "<<meta_file<<std::endl;
		return;
	}
	std::string str;
	std::string index;
	std::stringstream ss;
	MASS::Character* character = new MASS::Character();
	while(!ifs.eof())
	{
		str.clear();
		index.clear();
		ss.clear();

		std::getline(ifs,str);
		ss.str(str);
		ss>>index;
		if(!index.compare("use_muscle"))
		{	
			std::string str2;
			ss>>str2;
			if(!str2.compare("true"))
				this->SetUseMuscle(true);
			else
				this->SetUseMuscle(false);
		}
		else if(!index.compare("con_hz")){
			int hz;
			ss>>hz;
			this->SetControlHz(hz);
		}
		else if(!index.compare("sim_hz")){
			int hz;
			ss>>hz;
			this->SetSimulationHz(hz);
		}
		else if(!index.compare("sim_hz")){
			int hz;
			ss>>hz;
			this->SetSimulationHz(hz);
		}
		else if(!index.compare("skel_file")){
			std::string str2;
			ss>>str2;

			character->LoadSkeleton(std::string(MASS_ROOT_DIR)+str2,load_obj);
		}
		else if(!index.compare("muscle_file")){
			std::string str2;
			ss>>str2;
			if(this->GetUseMuscle())
				character->LoadMuscles(std::string(MASS_ROOT_DIR)+str2);
		}
		else if(!index.compare("bvh_file")){
			std::string str2,str3;

			ss>>str2>>str3;
			bool cyclic = false;
			if(!str3.compare("true"))
				cyclic = true;
			character->LoadBVH(std::string(MASS_ROOT_DIR)+str2,cyclic);
		}
		else if(!index.compare("reward_param")){
			double a,b,c,d;
			ss>>a>>b>>c>>d;
			this->SetRewardParameters(a,b,c,d);

		}


	}
	ifs.close();
	
	
	double kp = 300.0;
	character->SetPDParameters(kp,sqrt(2*kp));
	this->SetCharacter(character);
	this->SetGround(MASS::BuildFromFile(std::string(MASS_ROOT_DIR)+std::string("/data/ground.xml")));

	this->Initialize();
}
void
Environment::
Initialize()
{
	if(mCharacter->GetSkeleton()==nullptr){
		std::cout<<"Initialize character First"<<std::endl;
		exit(0);
	}
	if(mCharacter->GetSkeleton()->getRootBodyNode()->getParentJoint()->getType()=="FreeJoint")
		mRootJointDof = 6;
	else if(mCharacter->GetSkeleton()->getRootBodyNode()->getParentJoint()->getType()=="PlanarJoint")
		mRootJointDof = 3;	
	else
		mRootJointDof = 0;
	mNumActiveDof = mCharacter->GetSkeleton()->getNumDofs()-mRootJointDof;
	if(mUseMuscle)
	{
		int num_total_related_dofs = 0;
		for(auto m : mCharacter->GetMuscles()){
			m->Update();
			num_total_related_dofs += m->GetNumRelatedDofs();
		}
		mCurrentMuscleTuple.JtA = Eigen::VectorXd::Zero(num_total_related_dofs);
		mCurrentMuscleTuple.L = Eigen::MatrixXd::Zero(mNumActiveDof,mCharacter->GetMuscles().size());
		mCurrentMuscleTuple.b = Eigen::VectorXd::Zero(mNumActiveDof);
		mCurrentMuscleTuple.tau_des = Eigen::VectorXd::Zero(mNumActiveDof);
		mActivationLevels = Eigen::VectorXd::Zero(mCharacter->GetMuscles().size());
	}
	mWorld->setGravity(Eigen::Vector3d(0,-9.8,0.0));
	mWorld->setTimeStep(1.0/mSimulationHz);
	mWorld->getConstraintSolver()->setCollisionDetector(dart::collision::BulletCollisionDetector::create());
	mWorld->addSkeleton(mCharacter->GetSkeleton());
	mWorld->addSkeleton(mGround);
	mAction = Eigen::VectorXd::Zero(mNumActiveDof);
	
	Reset(false);
	mNumState = GetState().rows();
	
	for(std::size_t i = 0; i < mCharacter->GetSkeleton()->getNumRigidBodyNodes()-1 ; ++i)
        {    
        std::cout<<mCharacter->GetSkeleton()->getBodyNode(i)->getName()<<" "<<i<<std::endl;
        }	

	for(std::size_t i = 0; i < mCharacter->GetSkeleton()->getNumDofs() ; ++i)
        {    
        std::cout<<"getDof "<<mCharacter->GetSkeleton()->getDof(i)->getName()<<" "<<i<<std::endl;
        }
		
	reward_file.open ("1_reward_file.txt", std::ios::app);
	tibia_talus_aci.open ("2_tibia_talus_aci.txt", std::ios::app);
	tibia_yer_aci.open ("3_tibia_yer_aci.txt", std::ios::app);
   	femur_tibia_aci.open ("4_femur_tibia_aci.txt", std::ios::app);
	activationlevel_file.open ("5_activationlevel_file.csv", std::ios::app);
	head_neck_aci.open ("6_head_neck_aci.txt", std::ios::app);
	handR_acisal_hiz.open ("7_handR_acisal_hiz.txt", std::ios::app);
	handL_acisal_hiz.open ("8_handL_acisal_hiz.txt", std::ios::app);
	
	FemurL_x_desired_torque.open ("9_FemurL_x_desired_torque.txt", std::ios::app);
	FemurL_y_desired_torque.open ("10_FemurL_y_desired_torque.txt", std::ios::app);
	FemurL_z_desired_torque.open ("11_FemurL_z_desired_torque.txt", std::ios::app);
	TibiaL_desired_torque.open ("12_TibiaL_desired_torque.txt", std::ios::app);
	
	int count2 = 0;
	for(auto muscle : mCharacter->GetMuscles())
	{
		activationlevel_file<<mCharacter->GetMuscles()[count2++]->name<<",";
	}
	activationlevel_file<<std::endl;
	
	
}
void
Environment::
Reset(bool RSI)
{
	mWorld->reset();

	mCharacter->GetSkeleton()->clearConstraintImpulses();
	mCharacter->GetSkeleton()->clearInternalForces();
	mCharacter->GetSkeleton()->clearExternalForces();
	
	double t = 0.0;

	if(RSI)
		t = dart::math::random(0.0,mCharacter->GetBVH()->GetMaxTime()*0.9);
	mWorld->setTime(t);
	mCharacter->Reset();

	mAction.setZero();

	std::pair<Eigen::VectorXd,Eigen::VectorXd> pv = mCharacter->GetTargetPosAndVel(t,1.0/mControlHz);
	mTargetPositions = pv.first;
	mTargetVelocities = pv.second;

	mCharacter->GetSkeleton()->setPositions(mTargetPositions);
	mCharacter->GetSkeleton()->setVelocities(mTargetVelocities);
	mCharacter->GetSkeleton()->computeForwardKinematics(true,false,false);
	
	reward_file.close();
	activationlevel_file.close();
	tibia_yer_aci.close();
	femur_tibia_aci.close();
    	tibia_talus_aci.close();
	head_neck_aci.close();
	handR_acisal_hiz.close();
	handL_acisal_hiz.close();
	
	FemurL_x_desired_torque.close(); 
	FemurL_y_desired_torque.close();
	FemurL_z_desired_torque.close();
	TibiaL_desired_torque.close();
}


double x_tmp1 = 0;
double x_tmp2;
double x_tmp3 = 0;
double x_tmp4;
double x_tmp5;
double x_tmp6;
double x_tmp7;

void
Environment::
ProsthesisControl()
{

	reward_file<<GetReward()<<std::endl;  
	
	
}

void
Environment::
WriteActivation()
{
	int count1 = 0;
	for(auto muscle : mCharacter->GetMuscles())
	{
		activationlevel_file<<mActivationLevels[count1++]<<",";

	}
	activationlevel_file<<std::endl;
	//activationlevel_file<<mCharacter->GetMuscles()[0]->name<<",";
	
}
/////////////////////////////////////////////////////////////////////////////////////

//addExtForce(force, location, true, true);
//https://dartsim.github.io/tutorials_pendulum.html#lesson-2-set-spring-and-damping-properties-for-joints
/* void dart::dynamics::BodyNode::addExtForce 	( 	const Eigen::Vector3d &  	_force,
		const Eigen::Vector3d &  	_offset = Eigen::Vector3d::Zero(),
		bool  	_isForceLocal = false,
		bool  	_isOffsetLocal = true 
	) 	 */
/////////////////////////////////////////////////////////////////////////////////////
void
Environment::
disKuvvetUygulaOrta()
{
mCharacter->GetSkeleton()->getBodyNode("TibiaL")->addExtForce(
  Eigen::Vector3d(37.763, 0.21, 0), // linear force expressed in world coordinates
  Eigen::Vector3d(-0.07, 0.21, 0),  // offset from the origin of the body frame that you apply the linear force at. The offset is expressed in local coordinates
  true,
  true
);
}

void
Environment::
disKuvvetUygulaUst()
{
mCharacter->GetSkeleton()->getBodyNode("FemurL")->addExtForce(
  Eigen::Vector3d(-37.763, -0.1, 0), 
  Eigen::Vector3d(0.06, -0.1, 0),
  true,
  true  
);
}

void
Environment::
disKuvvetUygulaAlt()
{
mCharacter->GetSkeleton()->getBodyNode("TibiaL")->addExtForce(
  Eigen::Vector3d(-37.763, 0.1, 0), 
  Eigen::Vector3d(0.0575, 0.1, 0),
  true,
  true
);
}
/////////////////////////////////////////////////////////////////////////////////////
void
Environment::
disTorkUygulaOrta()
{
mCharacter->GetSkeleton()->getBodyNode("TibiaL")->addExtTorque(
  Eigen::Vector3d(1, 1, 1)  // torque (or angular force) expressed in world coordinates
);
}

void
Environment::
disTorkUygulaUst()
{
mCharacter->GetSkeleton()->getBodyNode("FemurL")->addExtTorque(
  Eigen::Vector3d(1, 1, 1)  
);
}

void
Environment::
disTorkUygulaAlt()
{
mCharacter->GetSkeleton()->getBodyNode("TibiaL")->addExtTorque(
  Eigen::Vector3d(1, 1, 1)  
);
}
/////////////////////////////////////////////////////////////////////////////////////

void
Environment::
disEtkiGosterOrta()
{
    ArrowShape::Properties arrow_properties;
    arrow_properties.mRadius = 0.008;
    mArrow = std::shared_ptr<ArrowShape>(new ArrowShape(
             //Eigen::Vector3d(0.0295, (0.1215*arrow_pow/0.1215) , -0.0695 ),
             //Eigen::Vector3d(0.0295, 0.0678 , -0.0695 ),
//vektorler arasindaki vektorel uzaklık uygula. ornek vektor1x=1 vector2x=3. x1 den x3e 2 birimlik. Fark 0sa o yonde vektor olusmaz.
  Eigen::Vector3d(-0.145, 0.21, 0), // linear force expressed in world coordinates (   //okun baslangic konumu +- degisikligi yon degistirir
  Eigen::Vector3d(-0.07, 0.21, 0),  // offset from the origin of the body frame that you apply the linear force at. The offset is expressed in local coordinates // okun ucunun konumu x+ saga,y+ yukariya,z+ sayfa disina 
             arrow_properties, dart::Color::Red(1.0)));  
      
    auto visualShapeNodes = mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getShapeNodesWith<VisualAspect>();
    
    visualShapeNodes[0]->getVisualAspect()->setColor(dart::Color::Green());   // "kutu rengi"
    
    mCharacter->GetSkeleton()->getBodyNode("TibiaL")->createShapeNodeWith<VisualAspect>(mArrow);
  //  mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getShapeNodesWith<VisualAspect>()[0]->getVisualAspect()->setColor(dart::Color::Blue());


//std::shared_ptr<ArrowShape> mArrow; >>>Environment:: Environment() altina
//arrow temizlemeyi unutma
// public: 	void disKuvvet(); environment.h a ekle
//Step() icine disKuvvet();
}

void
Environment::
disEtkiGosterUst()
{
    ArrowShape::Properties arrow_properties;
    arrow_properties.mRadius = 0.008;
    mArrow = std::shared_ptr<ArrowShape>(new ArrowShape(

  Eigen::Vector3d(0.135, -0.1, 0),
  Eigen::Vector3d(0.06, -0.1, 0),
             arrow_properties, dart::Color::Red(1.0)));  
      
    auto visualShapeNodes = mCharacter->GetSkeleton()->getBodyNode("FemurL")->getShapeNodesWith<VisualAspect>();
    
    visualShapeNodes[0]->getVisualAspect()->setColor(dart::Color::Green());   // "kutu rengi"
    
    mCharacter->GetSkeleton()->getBodyNode("FemurL")->createShapeNodeWith<VisualAspect>(mArrow);
}

void
Environment::
disEtkiGosterAlt()
{
    ArrowShape::Properties arrow_properties;
    arrow_properties.mRadius = 0.008;
    mArrow = std::shared_ptr<ArrowShape>(new ArrowShape(

  Eigen::Vector3d(0.1325, 0.1, 0),
  Eigen::Vector3d(0.0575, 0.1, 0),
             arrow_properties, dart::Color::Red(1.0)));  
      
    auto visualShapeNodes = mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getShapeNodesWith<VisualAspect>();
    
    visualShapeNodes[0]->getVisualAspect()->setColor(dart::Color::Green());   // "kutu rengi"
    
    mCharacter->GetSkeleton()->getBodyNode("TibiaL")->createShapeNodeWith<VisualAspect>(mArrow);
}
/////////////////////////////////////////////////////////////////////////////////////

void
Environment::
Step()
{	
	if(mUseMuscle)
	{
		int count = 0;
		for(auto muscle : mCharacter->GetMuscles())
		{
			muscle->activation = mActivationLevels[count++];
			muscle->Update();
			muscle->ApplyForceToBody();
		}
		if(mSimCount == mRandomSampleIndex)
		{
			auto& skel = mCharacter->GetSkeleton();
			auto& muscles = mCharacter->GetMuscles();

			int n = skel->getNumDofs();
			int m = muscles.size();
			Eigen::MatrixXd JtA = Eigen::MatrixXd::Zero(n,m);
			Eigen::VectorXd Jtp = Eigen::VectorXd::Zero(n);

			for(int i=0;i<muscles.size();i++)
			{
				auto muscle = muscles[i];
				// muscle->Update();
				Eigen::MatrixXd Jt = muscle->GetJacobianTranspose();
				auto Ap = muscle->GetForceJacobianAndPassive();

				JtA.block(0,i,n,1) = Jt*Ap.first;
				Jtp += Jt*Ap.second;
			}

			mCurrentMuscleTuple.JtA = GetMuscleTorques();
			mCurrentMuscleTuple.L = JtA.block(mRootJointDof,0,n-mRootJointDof,m);
			mCurrentMuscleTuple.b = Jtp.segment(mRootJointDof,n-mRootJointDof);
			mCurrentMuscleTuple.tau_des = mDesiredTorque.tail(mDesiredTorque.rows()-mRootJointDof);
			mMuscleTuples.push_back(mCurrentMuscleTuple);
		}
	}
	else
	{
		GetDesiredTorques();
		mCharacter->GetSkeleton()->setForces(mDesiredTorque);
	}


	
	//std::cout<<"getDof "<<mCharacter->GetSkeleton()->getDof(19)->getName()<<"  "<<mCharacter->GetSkeleton()->getDof(19)->getPosition()<<std::endl;
	//std::cout<<"getDof "<<mCharacter->GetSkeleton()->getDof(20)->getName()<<"  "<<mCharacter->GetSkeleton()->getDof(20)->getPosition()<<std::endl;
	//std::cout<<"getDof "<<mCharacter->GetSkeleton()->getDof(21)->getName()<<"  "<<mCharacter->GetSkeleton()->getDof(21)->getPosition()<<std::endl;
	
	WriteActivation();
		
	// TibiaL - Ground Angle  (Add current velocity to cumulative velocity)
	x_tmp1=matrixToEulerZXY(mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getWorldTransform().linear())(1);
	
	// FemurL - TibiaL Angle
	x_tmp2 = *(mCharacter->GetSkeleton()->getBodyNode("FemurL")->getChildJoint(0)->getPositions().data());
	
	// TibiaL Angular Velocity
	x_tmp3 = mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getAngularVelocity(Frame::World(),Frame::World())(0);
	
	// TibiaL - TalusL Angle
	x_tmp4 = *(mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getChildJoint(0)->getPositions().data());
	
	// Head - Neck Angle
	x_tmp5 = *(mCharacter->GetSkeleton()->getBodyNode("Neck")->getChildJoint(0)->getPositions().data());
	
	// handR Angular Velocity
	x_tmp6 = mCharacter->GetSkeleton()->getBodyNode("HandR")->getAngularVelocity(Frame::World(),Frame::World())(0);
	
	// handL Angular Velocity
	x_tmp7 = mCharacter->GetSkeleton()->getBodyNode("HandL")->getAngularVelocity(Frame::World(),Frame::World())(0);
	
	//std::cout<<"WORLDTRANSFORM"<<x_tmp1<<std::endl;

	tibia_yer_aci<<x_tmp1<<std::endl;
    	femur_tibia_aci<<x_tmp2<<std::endl;   
	tibia_talus_aci<<x_tmp4<<std::endl;  
	head_neck_aci<<x_tmp5<<std::endl;  
	handR_acisal_hiz<<x_tmp6<<std::endl;  
	handL_acisal_hiz<<x_tmp7<<std::endl;  
	
	mWorld->step();
	// Eigen::VectorXd p_des = mTargetPositions;
	// //p_des.tail(mAction.rows()) += mAction;
	// mCharacter->GetSkeleton()->setPositions(p_des);
	// mCharacter->GetSkeleton()->setVelocities(mTargetVelocities);
	// mCharacter->GetSkeleton()->computeForwardKinematics(true,false,false);
	// mWorld->setTime(mWorld->getTime()+mWorld->getTimeStep());
	//reward_file<<GetReward()<<std::endl;
	
	//disKuvvetUygulaOrta();
	//disKuvvetUygulaUst();
	//disKuvvetUygulaAlt();
	//disTorkUygulaOrta();
	//disTorkUygulaUst();
	//disTorkUygulaAlt();
	
	//std::cout<<"Force: "<<*(mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getChildJoint(0)->getExternalForces().data())<<std::endl;
	//std::cout<<"External Force: "<<mCharacter->GetSkeleton()->getBodyNode("TibiaL")->getExternalForceLocal()<<std::endl;
	
	FemurL_x_desired_torque<<mDesiredTorque[15]<<std::endl; 
	FemurL_y_desired_torque<<mDesiredTorque[16]<<std::endl; 
	FemurL_z_desired_torque<<mDesiredTorque[17]<<std::endl; 
	TibiaL_desired_torque<<mDesiredTorque[18]<<std::endl; 
/*
//std::cout<<mDesiredTorque[19];
FemurL_x 15
FemurL_y 16
FemurL_z 17
TibiaL 18

FemurL 6
TibiaL 7
*/
	
	mSimCount++;
}


Eigen::VectorXd
Environment::
GetDesiredTorques()
{
	Eigen::VectorXd p_des = mTargetPositions;
	p_des.tail(mTargetPositions.rows()-mRootJointDof) += mAction;
	mDesiredTorque = mCharacter->GetSPDForces(p_des);
	return mDesiredTorque.tail(mDesiredTorque.rows()-mRootJointDof);
}
Eigen::VectorXd
Environment::
GetMuscleTorques()
{
	int index = 0;
	mCurrentMuscleTuple.JtA.setZero();
	for(auto muscle : mCharacter->GetMuscles())
	{
		muscle->Update();
		Eigen::VectorXd JtA_i = muscle->GetRelatedJtA();
		mCurrentMuscleTuple.JtA.segment(index,JtA_i.rows()) = JtA_i;
		index += JtA_i.rows();
	}
	
	return mCurrentMuscleTuple.JtA;
}
double exp_of_squared(const Eigen::VectorXd& vec,double w)
{
	return exp(-w*vec.squaredNorm());
}
double exp_of_squared(const Eigen::Vector3d& vec,double w)
{
	return exp(-w*vec.squaredNorm());
}
double exp_of_squared(double val,double w)
{
	return exp(-w*val*val);
}


bool
Environment::
IsEndOfEpisode()
{
	bool isTerminal = false;
	
	Eigen::VectorXd p = mCharacter->GetSkeleton()->getPositions();
	Eigen::VectorXd v = mCharacter->GetSkeleton()->getVelocities();

	double root_y = mCharacter->GetSkeleton()->getBodyNode(0)->getTransform().translation()[1] - mGround->getRootBodyNode()->getCOM()[1];
	if(root_y<1.3)
		isTerminal =true;
	else if (dart::math::isNan(p) || dart::math::isNan(v))
		isTerminal =true;
	else if(mWorld->getTime()>10.0)
		isTerminal =true;
	
	return isTerminal;
}
Eigen::VectorXd 
Environment::
GetState()
{
	auto& skel = mCharacter->GetSkeleton();
	dart::dynamics::BodyNode* root = skel->getBodyNode(0);
	int num_body_nodes = skel->getNumBodyNodes() - 1;
	Eigen::VectorXd p,v;

	p.resize( (num_body_nodes-1)*3);
	v.resize((num_body_nodes)*3);

	for(int i = 1;i<num_body_nodes;i++)
	{
		p.segment<3>(3*(i-1)) = skel->getBodyNode(i)->getCOM(root);
		v.segment<3>(3*(i-1)) = skel->getBodyNode(i)->getCOMLinearVelocity();
	}
	
	v.tail<3>() = root->getCOMLinearVelocity();

	double t_phase = mCharacter->GetBVH()->GetMaxTime();
	double phi = std::fmod(mWorld->getTime(),t_phase)/t_phase;

	p *= 0.8;
	v *= 0.2;

	Eigen::VectorXd state(p.rows()+v.rows()+1);

	state<<p,v,phi;
	return state;
}
void 
Environment::
SetAction(const Eigen::VectorXd& a)
{
	mAction = a*0.1;

	double t = mWorld->getTime();

	std::pair<Eigen::VectorXd,Eigen::VectorXd> pv = mCharacter->GetTargetPosAndVel(t,1.0/mControlHz);
	mTargetPositions = pv.first;
	mTargetVelocities = pv.second;

	mSimCount = 0;
	mRandomSampleIndex = rand()%(mSimulationHz/mControlHz);
	mAverageActivationLevels.setZero();
}
double 
Environment::
GetReward()
{
	auto& skel = mCharacter->GetSkeleton();

	Eigen::VectorXd cur_pos = skel->getPositions();
	Eigen::VectorXd cur_vel = skel->getVelocities();
	
	//std::cout<<"cur_pos pre "<<cur_pos[20]<<std::endl;


	
	//std::cout<<"cur_pos post "<<cur_pos[20]<<std::endl;
	
	
	Eigen::VectorXd p_diff_all = skel->getPositionDifferences(mTargetPositions,cur_pos);
	Eigen::VectorXd v_diff_all = skel->getPositionDifferences(mTargetVelocities,cur_vel);
	
	p_diff_all[20]=0;
	p_diff_all[21]=0;
	p_diff_all[22]=0;
	p_diff_all[23]=0;
	v_diff_all[20]=0;
	v_diff_all[21]=0;
	v_diff_all[22]=0;
	v_diff_all[23]=0;

	Eigen::VectorXd p_diff = Eigen::VectorXd::Zero(skel->getNumDofs());
	Eigen::VectorXd v_diff = Eigen::VectorXd::Zero(skel->getNumDofs());

	const auto& bvh_map = mCharacter->GetBVH()->GetBVHMap();

	for(auto ss : bvh_map)
	{
		auto joint = mCharacter->GetSkeleton()->getBodyNode(ss.first)->getParentJoint();
		int idx = joint->getIndexInSkeleton(0);
		if(joint->getType()=="FreeJoint")
			continue;
		else if(joint->getType()=="RevoluteJoint")
			p_diff[idx] = p_diff_all[idx];
		else if(joint->getType()=="BallJoint")
			p_diff.segment<3>(idx) = p_diff_all.segment<3>(idx);
	}


	
	
	
	auto ees = mCharacter->GetEndEffectors();
	Eigen::VectorXd ee_diff(ees.size()*3);
	Eigen::VectorXd com_diff;

	
	
	
	for(int i =0;i<ees.size();i++)
		ee_diff.segment<3>(i*3) = ees[i]->getCOM();
	com_diff = skel->getCOM();

	skel->setPositions(mTargetPositions);
	skel->computeForwardKinematics(true,false,false);

	com_diff -= skel->getCOM();
	for(int i=0;i<ees.size();i++)
		ee_diff.segment<3>(i*3) -= ees[i]->getCOM()+com_diff;

	skel->setPositions(cur_pos);
	skel->computeForwardKinematics(true,false,false);
	
	
	Eigen::VectorXd p_diff_y(p_diff.head(20).size() + p_diff.tail(32).size());
	p_diff_y << p_diff.head(20), p_diff.tail(32);

	double r_q = exp_of_squared(p_diff_y,2.0);
	double r_v = exp_of_squared(v_diff,0.1);
	double r_ee = exp_of_squared(ee_diff,40.0);
	double r_com = exp_of_squared(com_diff,10.0);

	double r = r_ee*(w_q*r_q + w_v*r_v);
	
	
	//kuvvetin yerini gosterme amacli. Egitim sirasinda burasi kapalı olunca hizli oluyor. NN klasoru kaydet. Burayi ac. Calistir. NNleri yerine getir.
	//disEtkiGosterOrta();
	//disEtkiGosterUst();
	//disEtkiGosterAlt();	
	
	//std::cout<<"r_ee  "<<r_ee<<std::endl;
	return r;
}
