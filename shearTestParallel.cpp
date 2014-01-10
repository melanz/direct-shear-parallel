#include "common/common.h"
#include "common/input_output.h"

vector<char*> bodyTypes;

struct bodyLocation {
        double x;
        double y;
        double z;
        double e0;
        double e1;
        double e2;
        double e3;
};

int importBodyLocations(char* fileName, vector<bodyLocation>* bodyLocations)
{
        double pos_x=0;
        double pos_y=0;
        double pos_z=0;
        double rot_0=0;
        double rot_1=0;
        double rot_2=0;
        double rot_3=0;
		double vx = 0;
		double vy = 0;
		double vz = 0;
		int type = 0;
        std::string temp_data;
        ifstream ifile(fileName);

        bodyLocation bodyLoc;
        //for (int i = 0; i < numBodies; i++)
	while(getline(ifile,temp_data))
        {
                //getline(ifile,temp_data);
                for(int i=0; i<temp_data.size(); ++i){
                        if(temp_data[i]==','){temp_data[i]=' ';}
                }

                std::stringstream ss(temp_data);
                ss>>pos_x>>pos_y>>pos_z>>rot_0>>rot_1>>rot_2>>rot_3>>vx>>vy>>vz>>type;

                bodyLoc.x = pos_x;
                bodyLoc.y = pos_y;
                bodyLoc.z = pos_z;
                bodyLoc.e0 = rot_0;
                bodyLoc.e1 = rot_1;
                bodyLoc.e2 = rot_2;
                bodyLoc.e3 = rot_3;

                if(type==0) bodyLocations->push_back(bodyLoc);
        }

        return 0;
}

double getRandomNumber(double min, double max)
{
   // x is in [0,1[
   double x = rand()/static_cast<double>(RAND_MAX);

   // [0,1[ * (max - min) + min is in [min,max[
   double that = min + ( x * (max - min) );

   return that;
}

int generateRockObject(int rockIndex, double rockRadius)
{
        char filename[100];
        sprintf(filename, "../data/directShear/rocks/rock%d.obj", rockIndex);
        ChStreamOutAsciiFile rockFile(filename);

        rockFile << "# OBJ file created by ply_to_obj.c\n#\ng Object001";

        // add vertices
        //1
        ChVector<double> vertex = ChVector<>(-0.57735, -0.57735, 0.57735);
        ChVector<double> vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " <<vertex.z;
        //2
        vertex = ChVector<>(0.934172, 0.356822, 0);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //3
        vertex = ChVector<>(0.934172, -0.356822, 0);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //4
        vertex = ChVector<>(-0.934172, 0.356822, 0);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //5
        vertex = ChVector<>(-0.934172, -0.356822, 0);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //6
        vertex = ChVector<>(0, 0.934172, 0.356822);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //7
        vertex = ChVector<>(0, 0.934172, -0.356822);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //8
        vertex = ChVector<>(0.356822, 0, -0.934172);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //9
        vertex = ChVector<>(-0.356822, 0, -0.934172);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //10
        vertex = ChVector<>(0, -0.934172, -0.356822);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //11
        vertex = ChVector<>(0, -0.934172, 0.356822);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //12
        vertex = ChVector<>(0.356822, 0, 0.934172);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //13
        vertex = ChVector<>(-0.356822, 0, 0.934172);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //14
        vertex = ChVector<>(0.57735, 0.57735, -0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //15
        vertex = ChVector<>(0.57735, 0.57735, 0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //16
        vertex = ChVector<>(-0.57735, 0.57735, -0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //17
        vertex = ChVector<>(-0.57735, 0.57735, 0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //18
        vertex = ChVector<>(0.57735, -0.57735, -0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //19
        vertex = ChVector<>(0.57735, -0.57735, 0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;
        //20
        vertex = ChVector<>(-0.57735, -0.57735, -0.57735);
        vertPert = ChVector<>(getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3),getRandomNumber(-0.3,0.3));
        vertex = rockRadius*(vertex+vertPert);
        rockFile << "\nv " << vertex.x << " " << vertex.y << " " << vertex.z;

        // add faces
        rockFile << "\n\nf 19 3 2\nf 12 19 2\nf 15 12 2\nf 8 14 2\nf 18 8 2\nf 3 18 2\nf 20 5 4\nf 9 20 4\nf 16 9 4\nf 13 17 4\nf 1 13 4\nf 5 1 4\nf 7 16 4\nf 6 7 4\nf 17 6 4\nf 6 15 2\nf 7 6 2\nf 14 7 2\nf 10 18 3\nf 11 10 3\nf 19 11 3\nf 11 1 5\nf 10 11 5\nf 20 10 5\nf 20 9 8\nf 10 20 8\nf 18 10 8\nf 9 16 7\nf 8 9 7\nf 14 8 7\nf 12 15 6\nf 13 12 6\nf 17 13 6\nf 13 1 11\nf 12 13 11\nf 19 12 11\n";

        return 0;
}

ChSharedPtr<ChBody> createSphere(ChSystemParallel* mphysicalSystem, double radius, ChVector<> position, ChColor color, bool visualize, double mass, int collisionFamily, double friction)
{
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(friction);

        // create the body
        ChSharedPtr<ChBody> sphere = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
        //sphere->SetPos(position);
        //sphere->SetMass(mass);
        //mphysicalSystem->Add(sphere);

        // Define the collision shape
        //sphere->GetCollisionModel()->ClearModel();
        //sphere->GetCollisionModel()->AddSphere(radius);
        //sphere->GetCollisionModel()->BuildModel();
        //sphere->SetCollide(true);
	InitObject(sphere, mass, position, ChQuaternion<>(1, 0, 0, 0), material, true, false, collisionFamily, collisionFamily);
	AddCollisionGeometry(sphere, SPHERE, radius, ChVector<>(0,0,0), ChQuaternion<>(1,0,0,0));
	FinalizeObject(sphere, (ChSystemParallel *) mphysicalSystem);

        // Add visualization geometry
        //ChSharedPtr<ChSphereShape> sphereShape1(new ChSphereShape);
        //sphereShape1->GetSphereGeometry().rad = radius;
        //sphere->AddAsset(sphereShape1);

        // Add color
        //ChSharedPtr<ChVisualization> sphereColor(new ChVisualization);
        //sphereColor->SetColor(color);
        //sphere->AddAsset(sphereColor);

        return sphere;
}

void createParticlesFromFile(ChSystemParallel* mphysicalSystem, char* fileName, double particleRadius, double scalingFactor, double particleDensity, bool useSpheres, double friction, double visualize, ChColor color, double L, double H, double W, bool allOnes)
{
        vector<bodyLocation> bodyLocs;
        int numOtherBodies = importBodyLocations(fileName, &bodyLocs);

        ChSharedPtr<ChBody> particle;
        
        //ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
        //mmaterial->SetFriction(mu);

        // determine mass
        double volume = 4.0*CH_C_PI*particleRadius*particleRadius*particleRadius/3.0;
        double mass = particleDensity*volume;
        double Ix = .4*mass*particleRadius*particleRadius;
        if(allOnes)
        {
        	mass = 1.0;
        	Ix = 1.0;
        }

        int numRocks = 0;
        char filename[100];
        for(int i=numOtherBodies;i<bodyLocs.size();i++)
        {
                bodyLocation BL = bodyLocs[i];

                //if(BL.z<(.5*H-particleRadius)&&BL.z>(-.5*H+particleRadius))
                {

					if(!useSpheres)
					{
							//// USE ROCKS
							//sprintf(filename, "../data/directShear/rocks/rock%d.obj", numRocks);
							//particle = (ChBodySceneNode*)addChBodySceneNode_easyGenericMesh(
							//        &mphysicalSystem, msceneManager,
							//        mass,
							//        ChVector<>(BL.x,BL.y,BL.z),
							//        ChQuaternion<>(BL.e0,BL.e1,BL.e2,BL.e3),
							//        filename,
							//        false,        // not static
							//        false);        // true=convex; false=concave(do convex decomposition of concave mesh)
					}
					else
					{
							// USE SPHERES
							particle = createSphere(mphysicalSystem, particleRadius, ChVector<>(BL.x,BL.y,BL.z), color, mass, visualize, numRocks, friction);

							//particleEasy = (ChBodySceneNode*)addChBodySceneNode_easySphere(
							//        &mphysicalSystem, msceneManager,
							//        mass, // mass
							//        ChVector<>(BL.x,BL.y,BL.z),
							//        getRandomNumber(particleRadius,particleRadius), // radius
							//        20, // hslices, for rendering
							//        15); // vslices, for rendering
					}


                particle->SetInertiaXX(ChVector<>(Ix, Ix, Ix));
                particle->SetCollide(true);
                bodyTypes.push_back("particle");
                //particle->SetMaterialSurface(mmaterial);
                numRocks++;
                }
        }
}

int createParticles(ChSystemParallel* mphysicalSystem, double particleRadius, double scalingFactor, double particleDensity, ChVector<> size, bool useSpheres, double friction, bool visualize, ChColor color, bool allOnes)
{
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(0.4);

        double L = size.x;
        double H = size.y;
        double W = size.z;

        ChSharedPtr<ChBody> particle;
        //ChBodySceneNode* particleEasy;

//        ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
//        mmaterial->SetFriction(mu);

        // determine mass
        double volume = 4.0*CH_C_PI*particleRadius*particleRadius*particleRadius/3.0;
        double mass = particleDensity*volume;
        printf("Particle mass: %f\n",mass);
        double Ix = .4*mass*particleRadius*particleRadius;
        if(allOnes)
        {
        	mass = 1.0;
        	Ix = 1.0;
        }

        // make sure the particles are inside
        L = 0.6*L;
        W = 0.6*W;
        H = 15*H;

        int nLength = L/(2.4*particleRadius);
        int nWidth = W/(2.4*particleRadius);
        int nHeight = H/(2.4*particleRadius);
        printf("%d\n",nLength*nWidth*nHeight);

        //2.4*particleRadius*(H/(2.4*particleRadius))

        int numRocks = 0;
        char filename[100];
        for(int i=0;i<nLength;i++)//for(int i=0;i<1;i++)//
        {
        	for(int j=0;j<200;j++)//for(int j=0;j<nHeight;j++)//	
                {
                        for(int k=0;k<nWidth;k++)//for(int k=0;k<1;k++)//
                        {                        
                                if(!useSpheres)
                                {
                                        //// USE ROCKS
                                        //generateRockObject(numRocks,getRandomNumber(particleRadius,particleRadius));
                                        //sprintf(filename, "../data/directShear/rocks/rock%d.obj", numRocks);
                                        //particle = (ChBodySceneNode*)addChBodySceneNode_easyGenericMesh(
                                        //        &mphysicalSystem, msceneManager,
                                        //        mass,
                                        //        ChVector<>(2.4*particleRadius*i, 2.4*particleRadius*j, 2.4*particleRadius*k)+ChVector<>(-L/2+particleRadius, 1.4*particleRadius, -W/2+particleRadius)+ChVector<>(getRandomNumber(-.2*particleRadius,.2*particleRadius), getRandomNumber(-.2*particleRadius,.2*particleRadius), getRandomNumber(-.2*particleRadius,.2*particleRadius)), // pos,
                                        //        QUNIT,
                                        //        filename,
                                        //        false,        // not static
                                        //        false);        // true=convex; false=concave(do convex decomposition of concave mesh)
                                }
                                else
                                {
                                        // USE SPHERES
                                        particle = createSphere(mphysicalSystem, particleRadius,
                                                ChVector<>(2.4*particleRadius*i, 2.4*particleRadius*j, 2.4*particleRadius*k)+ChVector<>(-L/2+particleRadius, 2.4*particleRadius, -W/2+particleRadius),//+ChVector<>(getRandomNumber(-.2*particleRadius,.2*particleRadius), getRandomNumber(-.2*particleRadius,.2*particleRadius), getRandomNumber(-.2*particleRadius,.2*particleRadius)),
                                                //ChVector<>(0,1,0),//
                                                color, mass, visualize, numRocks, friction);
                                        //particle->GetCollisionModel()->SetFamily(5);
                                        //particle->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(5);

                                }

                                //particleEasy->GetBody()->SetInertiaXX(ChVector<>(Ix, Ix, Ix));
                                particle->SetInertiaXX(ChVector<>(Ix, Ix, Ix));
                                //particle->SetInertiaXX(particleEasy->GetBody()->GetInertiaXX());
                                //particle->SetCollide(true);
                                bodyTypes.push_back("particle");
                                //particle->SetMaterialSurface(mmaterial);
                                numRocks++;
				if(numRocks%10000==0) printf("Num rocks: %d\n",numRocks);
                        }
                }
        }

        return numRocks;
}

ChSharedPtr<ChBody> createShearPlate(ChSystemParallel* mphysicalSystem, ChVector<> size, double TH, ChVector<> position, ChQuaternion<> rotation, ChColor color, bool visualize, double scale, double friction)
{
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(friction);

        double L = size.x;
        double H = size.y;
        double W = size.z;
	double mass = 1;

        // create the body
	ChSharedPtr<ChBody> shearPlate = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

        // Define the collision shape
	InitObject(shearPlate, mass, position, rotation, material, true, false, -1, -1);
	AddCollisionGeometry(shearPlate, BOX, ChVector<>((L+2*TH)*scale*.5,H*scale*.5,TH*.5), ChVector<>(0,0,.5*W+.5*TH), rotation);
	AddCollisionGeometry(shearPlate, BOX, ChVector<>((L+2*TH)*scale*.5,H*scale*.5,TH*.5), ChVector<>(0,0,-.5*W-.5*TH), rotation);
	AddCollisionGeometry(shearPlate, BOX, ChVector<>(TH*.5,H*scale*.5,(W+2*TH)*scale*.5), ChVector<>(.5*L+.5*TH,0,0), rotation);
	AddCollisionGeometry(shearPlate, BOX, ChVector<>(TH*.5,H*scale*.5,(W+2*TH)*scale*.5), ChVector<>(-.5*L-.5*TH,0,0), rotation);
	FinalizeObject(shearPlate, (ChSystemParallel *) mphysicalSystem);

        return shearPlate;
}

ChSharedPtr<ChBody> createBox(ChSystemParallel* mphysicalSystem, ChVector<> size, ChVector<> position, ChQuaternion<> rotation, ChColor color, bool visualize, double friction)
{
	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(friction);

        double L = size.x;
        double H = size.y;
        double W = size.z;
	double mass = 1;

	// create the body
	ChSharedPtr<ChBody> box = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

        // Define the collision shape
	InitObject(box, mass, position, rotation, material, true, false, -1, -1);
	AddCollisionGeometry(box, BOX, ChVector<>(L*.5,H*.5,W*.5), ChVector<>(0,0,0), rotation);
	FinalizeObject(box, (ChSystemParallel *) mphysicalSystem);

        return box;
}

ChSharedPtr<ChBody> createShakerBox(ChSystemParallel* mphysicalSystem, ChVector<> size, double TH, double particleRadius, double friction)
{
        double L = size.x;
        double H = size.y;
        double W = size.z;

//        ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
//        mmaterial->SetFriction(friction); // Friction coefficient of steel

        ChSharedPtr<ChBody> ground = createBox(mphysicalSystem, ChVector<>((L+2*TH)*1.2, 2*TH, (W+2*TH)*1.2), ChVector<>(0,-.5*H+-.5*TH,0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false, friction);
       // ground->SetMaterialSurface(mmaterial);
        ground->SetBodyFixed(true);
        //ground->GetCollisionModel()->SetFamily(4);
        //ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
        bodyTypes.push_back("ground");

        // Create bottom
        ChSharedPtr<ChBody> bottom = createShearPlate(mphysicalSystem, ChVector<>(L,H*5,W), TH, ChVector<>(0, H*2-.5*H, 0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.3,0.6), false,1.2, friction);
        //bottom->SetMaterialSurface(mmaterial);
        bottom->SetBodyFixed(true);
        //bottom->GetCollisionModel()->SetFamily(4);
        //bottom->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
        bodyTypes.push_back("bottom");

        // Create cieling
        double height = 16*H+1.2*H;
        //ChSharedPtr<ChBody> cieling = createBox(mphysicalSystem, ChVector<>(L*.98, 2*H, W*.98), ChVector<>(0,height,0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false);
        //ChSharedPtr<ChBody> cieling = createBox(mphysicalSystem, ChVector<>(L, 2*H, W), ChVector<>(0,height,0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false);
        //cieling->SetMaterialSurface(mmaterial);
        //bodyTypes.push_back("cieling");
/*
        // create the translational joint between the top shear box and weight load
        ChVector<> groundCM = ground->GetPos();
        ChLinkLockLock* lock = new ChLinkLockLock();
        lock->Initialize(ground, bottom, ChCoordsys<>(groundCM, chrono::Q_from_AngAxis(-CH_C_PI/2,VECT_Y)) );
        mphysicalSystem.AddLink(lock);

        // apply motion
        ChFunction_Sine* motionFunc = new ChFunction_Sine(0,100,W*0.01);
        lock->SetMotion_Z(motionFunc);
*/
        return bottom;//cieling;
}

class ShearBox {
public:
        // THE DATA
        bool visualize;
        double L;
        double H;
        double W;
        double TH;
        double desiredVelocity;
        double normalPressure;

        // bodies
        ChSharedPtr<ChBody> ground;
        ChSharedPtr<ChBody> sideWall1;
        ChSharedPtr<ChBody> sideWall2;
        ChSharedPtr<ChBody> top;
        ChSharedPtr<ChBody> bottom;
        ChSharedPtr<ChBody> cieling;
        ChLinkLockLock* translational;

        // THE FUNCTIONS

        // Build and initialize the car, creating all bodies corresponding to
        // the various parts and adding them to the physical system - also creating
        // and adding constraints to the system.
        ShearBox(ChSystemParallel* mphysicalSystem,        ///< the chrono::engine physical system
                bool visualize, double L, double H, double W, double TH, double desiredVelocity, double normalPressure, double friction)
        {
                this->visualize = visualize;
                this->L = L;
                this->H = H;
                this->W = W;
                this->TH = TH;
                this->desiredVelocity = desiredVelocity;
                this->normalPressure = normalPressure;

//                ChSharedPtr<ChMaterialSurface> mmaterial(new ChMaterialSurface);
//                mmaterial->SetFriction(mu); // Friction coefficient of steel

                // Create ground
                ground = createBox(mphysicalSystem, ChVector<>((L+2*TH)*3,2*TH,(W+2*TH)*3), ChVector<>(0,-.5*H+-.5*TH,0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false, friction);
                //ground->SetMaterialSurface(mmaterial);
                ground->SetBodyFixed(true);
                //ground->GetCollisionModel()->SetFamily(4);
                //ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
                bodyTypes.push_back("ground");

                // Create side wall
                sideWall1 = createBox(mphysicalSystem, ChVector<>((L+2*TH)*3,H*3,TH), ChVector<>(0,0,.5*W+.5*TH), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false, friction);
                //sideWall1->SetMaterialSurface(mmaterial);
                sideWall1->SetBodyFixed(true);
                //ground->GetCollisionModel()->SetFamily(4);
                //ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
                bodyTypes.push_back("sideWall1");

                // Create side wall
                sideWall2 = createBox(mphysicalSystem, ChVector<>((L+2*TH)*3,H*3,TH), ChVector<>(0,0,-.5*W-.5*TH), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.6,0.6), false, friction);
                //sideWall2->SetMaterialSurface(mmaterial);
                sideWall2->SetBodyFixed(true);
                //ground->GetCollisionModel()->SetFamily(4);
                //ground->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
                bodyTypes.push_back("sideWall2");

                // Create bottom
                double shearPlateScale = .7;
                bottom = createShearPlate(mphysicalSystem, ChVector<>(L,shearPlateScale*H,W), 5*TH, ChVector<>(0, -.5*H+-.5*TH+shearPlateScale*.5*H+TH, 0), ChQuaternion<>(1,0,0,0), ChColor(0.6,0.3,0.6), false,1,friction);
                //bottom->SetMaterialSurface(mmaterial);
                //bottom->GetCollisionModel()->SetFamily(4);
                //bottom->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
                bottom->SetBodyFixed(true);
                bodyTypes.push_back("bottom");

                // Create top
                top = createShearPlate(mphysicalSystem, ChVector<>(L,shearPlateScale*H,W), 5*TH, ChVector<>(0, shearPlateScale*H+-.5*H+-.5*TH+shearPlateScale*.5*H+TH, 0), ChQuaternion<>(1,0,0,0), ChColor(0.3,0.3,0.6), false,1,friction);
                //top->SetMaterialSurface(mmaterial);
                //top->GetCollisionModel()->SetFamily(4);
               	//top->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
		//top->SetBodyFixed(true);
                bodyTypes.push_back("top");
/*
                // Create ground
                cieling = createBox(mphysicalSystem, ChVector<>(L, 2*H, W), ChVector<>(0,6.36,0), ChQuaternion<>(1,0,0,0), ChColor(0.3,0.5,0.3), false);
                cieling->SetMaterialSurface(mmaterial);
                //cieling->GetCollisionModel()->SetFamily(4);
                //cieling->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(4);
                bodyTypes.push_back("cieling");

                // apply force to cieling (SHOULD BE DONE BY WEIGHT)
                cieling->Empty_forces_accumulators();
                cieling->Accumulate_force(ChVector<>(0, -normalPressure*L*W, 0), cieling->GetPos(), false);

                // create the translational joint between the truss and weight load
                ChVector<> groundCM = ground->GetPos();
                ChSharedPtr<ChLinkLockPrismatic> translationalVert(new ChLinkLockPrismatic);
                translationalVert->Initialize(cieling, top,
                        ChCoordsys<>(groundCM, chrono::Q_from_AngAxis(CH_C_PI/2,VECT_X)) );
                mphysicalSystem->AddLink(translationalVert);
*/

                // create the translational joint between the top shear box and weight load
                ChVector<> topCM = top->GetPos();
                translational = new ChLinkLockLock();
                translational->Initialize(ground, top,
                        ChCoordsys<>(topCM, chrono::Q_from_AngAxis(-CH_C_PI/2,VECT_Y)) );
                mphysicalSystem->AddLink(translational);

                // apply motion
                ChFunction_Ramp* motionFunc3 = new ChFunction_Ramp(0,desiredVelocity);
                translational->SetMotion_Z(motionFunc3);

        }
};

template<class T>
void RunTimeStep(T* mSys, const int frame) {
	//mSys->DoStepDynamics(.00001);
}

int main(int argc, char* argv[]) 
{
	// command line arguments
	bool settleBodies = false;
	settleBodies = atoi(argv[1]); //if(argc==2) settleBodies = atoi(argv[1]);
	double contactRecoverySpeed = 100;
	contactRecoverySpeed = atoi(argv[2]); //if(argc==3) contactRecoverySpeed = atoi(argv[2]);
	int max_iteration = 30;
	max_iteration = atoi(argv[3]); //if(argc==4) max_iteration = atoi(argv[3]);
	real timestep = 0.0001;		// step size
	timestep = atof(argv[4]); //if(argc==5) timestep = atof(argv[4]);
	string data_folder = "./DataAllOnes/";
	data_folder = argv[5]; //if(argc==6) data_folder = argv[5];
	bool allOnes = atoi(argv[6]);

	//stringstream ss_dataFolder;
	//ss_dataFolder << "./data_test/data_" << contactRecoverySpeed << "_" << max_iteration << "_" << timestep;
	//string data_folder = ss_dataFolder.str();
	cout << "DATA FOLDER: " << data_folder << endl;

	bool visualize = false;
	int threads = 16;
	int config = 0;
	//real gravity = -9.81;			// acceleration due to gravity
	//real time_to_run = 1;			// length of simulation
	real current_time = 0;

	double tolerance = 1e-3;

	//=========================================================================================================
	// Create system
	//=========================================================================================================
	ChSystemParallel* mphysicalSystem = new ChSystemParallel();

	//=========================================================================================================
	// Populate the system with bodies/constraints/forces/etc.
	//=========================================================================================================

	bool useSpheres = true;
	double scalingFactor = 1000; // scaling for distance
	double scalingMASS = 10000000; // scaling for mass
	double L = .06*scalingFactor;
	double H = .03*scalingFactor;
	double W = .06*scalingFactor;//.06*scalingFactor;
	double TH = .0032*scalingFactor; // If this is changed I think the particle input might not be in the right place...
	//bool visualize = false;
	double desiredVelocity = .66e-3*scalingFactor;
	double particleRadius = .0004*scalingFactor;//.006*scalingFactor;//
	bool importParticles = true;
	double lengthToRun = .006*scalingFactor;
	double time_to_run = 5;//lengthToRun/desiredVelocity;
	int num_steps = time_to_run / timestep;
	double particleDensity = 2800*scalingMASS/scalingFactor/scalingFactor/scalingFactor;
	ShearBox* shearBox;
	ChSharedPtr<ChBody> cielingSettled;
	double muWalls = .2;
	double muParticles = .86;
	double gravity = -9.81*scalingFactor;

	// apply normal force to cieling
	// 16,888.1 Pa
	// 44,127.0 Pa
	// 71,365.9 Pa
	double normalPressure = 16888.1*scalingMASS/scalingFactor;
	//if(!settleBodies) gravity = 0;
	mphysicalSystem->Set_G_acc(ChVector<>(0,gravity,0));

	int numRocksCreated = 0;
	shearBox = new ShearBox(mphysicalSystem, visualize, L, H, W, TH, desiredVelocity, normalPressure, muWalls);
	if(!settleBodies)
	{
		createParticlesFromFile(mphysicalSystem, "posStart.txt", particleRadius, scalingFactor, particleDensity, useSpheres, muParticles, visualize, ChColor(0.6,0.6,0.6), L, H, W, allOnes);
		//shearBox = new ShearBox(mphysicalSystem, visualize, L, H, W, TH, desiredVelocity, normalPressure, muWalls);
	}
	else
	{
		ChVector<> size = chrono::ChVector<>(L,H,W);
		//cielingSettled  = createShakerBox(mphysicalSystem, size, TH, particleRadius, muWalls);
		shearBox->top->SetBodyFixed(true);
		numRocksCreated = createParticles(mphysicalSystem, particleRadius, scalingFactor, particleDensity, size, useSpheres, muParticles, visualize, ChColor(0.6,0.6,0.6), allOnes);
	}

	//=========================================================================================================
	// Edit system settings
	//=========================================================================================================

	mphysicalSystem->SetIntegrationType(ChSystem::INT_ANITESCU);
	mphysicalSystem->SetParallelThreadNumber(threads);
	mphysicalSystem->SetMaxiter(max_iteration);
	mphysicalSystem->SetIterLCPmaxItersSpeed(max_iteration);
	mphysicalSystem->SetTol(tolerance);
	mphysicalSystem->SetTolSpeeds(tolerance);
	//system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	mphysicalSystem->SetStep(timestep);

	((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->SetMaxIteration(max_iteration);
	((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->SetCompliance(0, 0, 0);
	((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->SetContactRecoverySpeed(300);
	((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->SetSolverType(ACCELERATED_PROJECTED_GRADIENT_DESCENT);

	((ChCollisionSystemParallel*) (mphysicalSystem->GetCollisionSystem()))->SetCollisionEnvelope(particleRadius * .05);
	((ChCollisionSystemParallel*) (mphysicalSystem->GetCollisionSystem()))->setBinsPerAxis(I3(25, 25, 25));
	((ChCollisionSystemParallel*) (mphysicalSystem->GetCollisionSystem()))->setBodyPerBin(100, 50);

	((ChSystemParallel*) mphysicalSystem)->SetAABB(R3(-L, -H, -W), R3(L, H, W));

	omp_set_num_threads(threads);

	//=========================================================================================================
	// Enter the time loop and render the simulation
	//=========================================================================================================
	if (visualize) {
		ChOpenGLManager * window_manager = new ChOpenGLManager();
		ChOpenGL openGLView(window_manager, mphysicalSystem, 800, 600, 0, 0, "Test_Solvers");
		openGLView.render_camera->camera_pos = Vector(0, 5, -20);
		openGLView.render_camera->look_at = Vector(0, 0, 0);
		openGLView.SetCustomCallback(RunTimeStep);
		openGLView.StartSpinning(window_manager);
		window_manager->CallGlutMainLoop();
	}

	//=========================================================================================================
	// If you choose not to visualize, output data for post-processing
	//=========================================================================================================

    // initialize the shear stress vs. displacement output file

	stringstream ss_m;
	ss_m << data_folder << "/" << "timing.txt";
	string timing_file_name = ss_m.str();
	ofstream ofile(timing_file_name.c_str());
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		//ofile.open();
		ofile << current_time << ", " << shearBox->top->GetPos().x << ", " << shearBox->top->GetPos().y << ", " << shearBox->top->GetPos().z << ", " << shearBox->translational->Get_react_force().x << ", " << shearBox->translational->Get_react_force().y << ", " << shearBox->translational->Get_react_force().z << ", " << ((ChLcpSolverParallel *) (mphysicalSystem->GetLcpSolverSpeed()))->GetResidual() << ", " << endl;

		cout << "step " << i;
		cout << " Residual: " << ((ChLcpSolverParallel *) (mphysicalSystem->GetLcpSolverSpeed()))->GetResidual();
		cout << " ITER: " << ((ChLcpSolverParallel *) (mphysicalSystem->GetLcpSolverSpeed()))->GetTotalIterations();
		cout << " OUTPUT STEP: Time= " << current_time << " bodies= " << mphysicalSystem->GetNbodies() << " contacts= " << mphysicalSystem->GetNcontacts() << " step time=" << mphysicalSystem->GetTimerStep()
				<< " lcp time=" << mphysicalSystem->GetTimerLcp() << " CDbroad time=" << mphysicalSystem->GetTimerCollisionBroad() << " CDnarrow time=" << mphysicalSystem->GetTimerCollisionNarrow() << " Iterations="
				<< ((ChLcpSolverParallel*) (mphysicalSystem->GetLcpSolverSpeed()))->GetTotalIterations() << "\n";
		//TimingFile(system_gpu, timing_file_name, current_time);
		mphysicalSystem->DoStepDynamics(timestep);
		RunTimeStep(mphysicalSystem, i);
		int save_every = 1.0 / timestep / 100.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(mphysicalSystem, ss.str());
			//output.ExportData(ss.str());
			file++;
		}
		current_time += timestep;
	}
	ofile.close();
	stringstream ss;
	return 0;
}

