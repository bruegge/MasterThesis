#include "Physics.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

void CPhysics::SetL1(double l1)
{
	m_dL1 = l1;
}
void CPhysics::SetL2(double l2)
{
	m_dL2 = l2;
}
void CPhysics::SetM1(double m1)
{
	m_dM1 = m1;
}
void CPhysics::SetM2(double m2)
{
	m_dM2 = m2;
}
void CPhysics::SetK(double k)
{
	m_dK = k;
}



double CPhysics::GetL1()
{
	return m_dL1;
}
double CPhysics::GetL2()
{
	return m_dL2;
}
double CPhysics::GetM1()
{
	return m_dM1;
}
double CPhysics::GetM2()
{
	return m_dM2;
}
double CPhysics::GetK()
{
	return m_dK;
}

void CPhysics::RecordState()
{
	SState state;

	switch (m_eDGL)
	{
	case CPhysics::PendulumAngular:
		state.vecVariables.resize(1);
		break;
	case CPhysics::PendulumSpring:
		state.vecVariables.resize(6);
		break;
	case CPhysics::PendulumPosition:
		state.vecVariables.resize(8);
		break;
	case CPhysics::LorenzAttractor:
		state.vecVariables.resize(3);
		break;
	case CPhysics::SinglePoint:
		state.vecVariables.resize(4);
		break;
	case CPhysics::EasyExample:
		state.vecVariables.resize(1);
		break;
	case CPhysics::DoublePendulum:
		state.vecVariables.resize(4);
		break;
	default:
		break;
	}

	if (m_eDGL == EDGL::PendulumPosition)
	{
		state.vecVariables[0] = m_sParticle.vecVariables[0];
		state.vecVariables[1] = m_sParticle.vecVariables[1];
		state.vecVariables[2] = m_sParticle.vecVariables[2];
		state.vecVariables[3] = m_sParticle.vecVariables[3];
		state.vecVariables[4] = m_sParticle.vecVariables[4];
		state.vecVariables[5] = m_sParticle.vecVariables[5];
		state.vecVariables[6] = m_sParticle.vecVariables[6];
		state.vecVariables[7] = m_sParticle.vecVariables[7];
		m_vecRecordedData.push_back(state);
	}
	else if (m_eDGL == EDGL::SinglePoint)
	{
		state.vecVariables[0] = m_sParticle.vecVariables[0];
		state.vecVariables[1] = m_sParticle.vecVariables[1];
		state.vecVariables[2] = m_sParticle.vecVariables[2];
		state.vecVariables[3] = m_sParticle.vecVariables[3];
		m_vecRecordedData.push_back(state);
	}
	else
	{
		for (unsigned int i = 0; i < state.vecVariables.size()- m_vAcceleration.size(); ++i)
		{
			state.vecVariables[i] = m_sParticle.vecVariables[i];
		}
		for (unsigned int i = 0; i < m_vAcceleration.size(); ++i)
		{
			state.vecVariables[m_sParticle.vecVariables.size() + i] = m_vAcceleration[i];
		}
		m_vecRecordedData.push_back(state);
		
	}
}

std::vector<long double> CPhysics::CalculateNewDerivation(std::vector<long double> particle)
{
	//Particle 2D
	double g = 9.807; //gravitations konstante  "Pendulum"
	double l = m_dL1; //länge des pendels
	double l2 = m_dL2; //länge des 2nd pendels
	double m1 = m_dM1; //mass 1
	double m2 = m_dM2; //mass 2
	double sigma = 10.0;
	double beta = 8.0 / 3.0;
	double rho = 28.0;
	double constant = 5.0;
	double k = m_dK; //springConst
	
	std::vector<long double> derivativeResult;
	derivativeResult.resize(particle.size());
	//PendulumSpring
	glm::vec2 position = glm::vec2(particle[0], particle[1]);
	glm::vec2 force = -glm::normalize(position) * static_cast<float>(k* (glm::length(position) - l)) + glm::vec2(0, -g);

	switch (m_eDGL)
	{
	case CPhysics::PendulumAngular:

		//first value: Theta
		//second value: U = (velocity)
		derivativeResult[0] = particle[1];
		derivativeResult[1] = -g / l * sinl(particle[0]);


		break;

	case CPhysics::PendulumSpring:

		derivativeResult[0] = particle[2];
		derivativeResult[1] = particle[3];


		derivativeResult[2] = force.x;
		derivativeResult[3] = force.y;
		m_vAcceleration[0] = force.x;
		m_vAcceleration[1] = force.y;
		break;
	case CPhysics::LorenzAttractor:

		derivativeResult[0] = sigma * (particle[1] - particle[0]);
		derivativeResult[1] = particle[0] * (rho - particle[2]) - particle[1];
		derivativeResult[2] = particle[0] * particle[1] - beta * particle[2];

		break;
	case CPhysics::SinglePoint:
		derivativeResult[0] = particle[2];
		derivativeResult[1] = particle[3];
		derivativeResult[2] = 0;
		derivativeResult[3] = -g;


		break;
	case CPhysics::EasyExample:
		derivativeResult[0] = 3;
	
		break;
	
	default:
		break;
	}

	if (m_eDGL == CPhysics::DoublePendulum)
	{
		double pT1 = particle[2];
		double pT2 = particle[3];
		double T1 = particle[0];
		double T2 = particle[1];

		double C1 = (pT1*pT2*sin(T1 - T2)) / (l*l2*(m1 + m2*sin(T1 - T2)*sin(T1 - T2)));
		double C2 = (l2*l2*m2*pT1*pT1 + l*l*(m1 + m2)*pT2*pT2 - l*l2*m2*pT1*pT2*cos(T1 - T2)) * sin(2 * (T1 - T2)) / (2 * l*l*l2*l2*pow(m1 + m2*sin(T1 - T2) * sin(T1 - T2), 2));

		derivativeResult[0] = (l2*pT1 - l*pT2 * cos(T1 - T2)) / (l*l*l2*(m1 + m2*sin(T1 - T2)*sin(T1 - T2))); //angle theta1
		derivativeResult[1] = (l*(m1 + m2)*pT2 - l2*m2*pT1*cos(T1 - T2)) / (l*l2*l2*m2*(m1 + m2*sin(T1 - T2)*sin(T1 - T2))); //angle theta2
		derivativeResult[2] = -(m1 + m2)*g*l*sin(T1) - C1 + C2; //momentum angle1
		derivativeResult[3] = -m2*g*l2*sin(T2) + C1 - C2; //momentum angle2
	}
	return derivativeResult;
}

void CPhysics::SetInitialValues(float h)
{
	m_vecRecordedData.clear();
	
	srand(time(NULL));
	std::stringstream ss;
	ss << "../../Generated trajectories/";
	
	switch (m_eDGL)
	{
	case CPhysics::PendulumAngular:
		ss << "pendulumAngular_L"<<m_dL1<<"_h"<<h<<"_";
		break;
	case CPhysics::PendulumPosition:
		ss << "pendulumPosition_L"<<m_dL1<<"_h"<<h << "_";
		break;
	case CPhysics::PendulumSpring:
		ss << "pendulumSpring_L"<<m_dL1<<"_h"<< h << "_";
		break;
	case CPhysics::LorenzAttractor:
		ss << "LorenzAttractor_h" << h << "_";
		break;
	case CPhysics::SinglePoint:
		ss << "SinglePoint_";
		break;
	case CPhysics::EasyExample:
		ss << "EasyExample_";
		break;
	case CPhysics::DoublePendulum:
		ss << "DoublePendulum_L1_" << m_dL1 << "_L2_"<<m_dL2<<"_m1_"<<m_dM1<<"_m2_"<<m_dM2<<"_h" << h<<"_";
		break;
	default:
		break;
	}

	for (unsigned int i = 0; i < 6; ++i)
	{
		ss << rand() % 10;
	}
	ss << ".txt";
	m_strFileName = ss.str();

	
	switch (m_eDGL)
	{
	case CPhysics::PendulumAngular:
		m_sParticle.vecVariables.resize(2);
	
		m_sParticle.vecVariables[0] = 3.1415926 / 180 * (85+90);  //Theta
		m_sParticle.vecVariables[1] = 0; //U = (velocity)

		break;
	case CPhysics::PendulumPosition:
		m_sParticle.vecVariables.resize(8);
	
		m_sParticle.vecVariables[0] = cos(3.1415926/180 * 85)*m_dL1;  //xPosition
		m_sParticle.vecVariables[1] = sin(3.1415926/180 * 85)*m_dL1;  //yPosition
		m_sParticle.vecVariables[2] = 0;  //xVelocity
		m_sParticle.vecVariables[3] = 0;  //yVelocity
		m_sParticle.vecVariables[4] = 0;  //xAcceleration
		m_sParticle.vecVariables[5] = 0;  //yAcceleration
		m_sParticle.vecVariables[6] = 0;  //correctionXPosition
		m_sParticle.vecVariables[7] = 0;  //correctionYPosition

		break;

	case CPhysics::PendulumSpring:
		m_sParticle.vecVariables.resize(4);
		m_vAcceleration.resize(2);
		m_sParticle.vecVariables[0] = cos(3.1415926 / 180 * 85) * m_dL1;  //xPosition
		m_sParticle.vecVariables[1] = sin(3.1415926 / 180 * 85) * m_dL1;  //yPosition
		m_sParticle.vecVariables[2] = 0;  //xVelocity
		m_sParticle.vecVariables[3] = 0;  //yVelocity

		break;
	case CPhysics::LorenzAttractor:
		m_sParticle.vecVariables.resize(3);
		
		m_sParticle.vecVariables[0] = 1.01;
		m_sParticle.vecVariables[1] = 1;
		m_sParticle.vecVariables[2] = 1;
		break;
	case CPhysics::SinglePoint:
		m_sParticle.vecVariables.resize(4);
		
		m_sParticle.vecVariables[0] = 0; //PosX
		m_sParticle.vecVariables[1] = 0; //PosY
		m_sParticle.vecVariables[2] = 5; //VeloX
		m_sParticle.vecVariables[3] = 20; //VeloY
		break;
	case CPhysics::EasyExample:
		m_sParticle.vecVariables.resize(1);
		
		m_sParticle.vecVariables[0] = 70; //PosX
		break;
	case CPhysics::DoublePendulum:
		m_sParticle.vecVariables.resize(4);
		
		m_sParticle.vecVariables[0] = 3.0; 
		m_sParticle.vecVariables[1] = 3.0 + 3.1415926/180.0;
		m_sParticle.vecVariables[2] = 0.0;
		m_sParticle.vecVariables[3] = 0.0;
	
		break;
	default:
		break;
	}
}

void CPhysics::SetCountRecordedStates(int nCount)
{
	m_nCountRecordedStates = nCount;
}

void CPhysics::SaveTrackedValues()
{
	std::ofstream myfile;
	myfile.precision(17);
	myfile.open(m_strFileName);
	for (unsigned int i = 0; i < m_vecRecordedData.size(); ++i)
	{
		for (unsigned int j = 0; j < m_vecRecordedData[i].vecVariables.size(); ++j)
		{
			myfile << m_vecRecordedData[i].vecVariables[j] << " ";
		}

		myfile << "\n";
	}
	myfile.close();

}

void CPhysics::SaveCurrentParticlePosition()
{
	//Pendulum
	float l = m_dL1; 
	float l2 = m_dL2;
	glm::mat4 rotationMatrix;
	glm::mat4 rotationMatrix2; 
	glm::mat4 translationMatrix;

	switch (m_eDGL)
	{
	case CPhysics::PendulumAngular:
		rotationMatrix = glm::rotate(rotationMatrix, static_cast<float>(m_sParticle.vecVariables[0]), glm::vec3(0, 0, 1));
		m_vParticlePosition = rotationMatrix * glm::vec4(0, -l, 0, 1);
		break;
	case CPhysics::PendulumPosition:
		m_vParticlePosition = glm::vec4(m_sParticle.vecVariables[0], m_sParticle.vecVariables[1], 0, 1);
		break;
	case CPhysics::PendulumSpring:
		m_vParticlePosition = glm::vec4(m_sParticle.vecVariables[0], m_sParticle.vecVariables[1], 0, 1);
		break;
	case CPhysics::LorenzAttractor:
		m_vParticlePosition = glm::vec4(m_sParticle.vecVariables[0], m_sParticle.vecVariables[1], 0, 1);
		break;
	case CPhysics::SinglePoint:
		m_vParticlePosition = glm::vec4(m_sParticle.vecVariables[0], m_sParticle.vecVariables[1], 0, 1);
		break;
	case CPhysics::EasyExample:
		m_vParticlePosition = glm::vec4(m_sParticle.vecVariables[0], 0, 0, 1);
		break;
	case CPhysics::DoublePendulum:
		rotationMatrix = glm::rotate(rotationMatrix, static_cast<float>(m_sParticle.vecVariables[0]), glm::vec3(0, 0, 1));
		m_vParticlePosition = rotationMatrix * glm::vec4(0, -l, 0, 1);
		rotationMatrix2 = glm::rotate(rotationMatrix2, static_cast<float>(m_sParticle.vecVariables[1]), glm::vec3(0, 0, 1));
		translationMatrix = glm::translate(translationMatrix, glm::vec3(m_vParticlePosition));
		m_vParticlePosition2 = translationMatrix * rotationMatrix2*glm::vec4(0, -l2, 0, 1);
		break;
	default:
		break;
	}

}

void CPhysics::SetDGL(const char* pDGL)
{
	if (pDGL == "PendulumAngular")
	{
		m_eDGL = EDGL::PendulumAngular;
	}
	else if (pDGL == "PendulumPosition")
	{
		m_eDGL = EDGL::PendulumPosition;
	}
	else if (pDGL == "PendulumSpring")
	{
		m_eDGL = EDGL::PendulumSpring;
	}
	else if (pDGL == "LorenzAttractor")
	{
		m_eDGL = EDGL::LorenzAttractor;
	}
	else if (pDGL == "SinglePoint")
	{
		m_eDGL = EDGL::SinglePoint;
	}
	else if (pDGL == "EasyExample")
	{
		m_eDGL = EDGL::EasyExample;
	}
	else if (pDGL == "DoublePendulum")
	{
		m_eDGL = EDGL::DoublePendulum;
	}
}


CPhysics::CPhysics()
{
	m_pModel = new CModel("../Models/Sphere.obj");

	glBindVertexArray(0);
	glGenVertexArrays(1, &m_nLineVAO);
	glBindVertexArray(m_nLineVAO);
	glGenBuffers(1, &m_nLineVBO);

	// load data into vertex buffers
	glBindBuffer(GL_ARRAY_BUFFER, m_nLineVBO);

	
	GLfloat vertices[] =
	{
		0, 0, 0, 0,
		1, 1, 1, 1,
		2, 2, 2, 2
	};


	glBufferData(GL_ARRAY_BUFFER, 3 * sizeof(glm::vec4), vertices, GL_STATIC_DRAW);


	glBindVertexArray(0);
	
}

CPhysics::~CPhysics()
{
	delete m_pModel;
}

void CPhysics::DrawParticle(CShader* pShader, CShader* pShaderLine, CCamera* pCamera)
{
	SaveCurrentParticlePosition();
	pShader->Bind();

	glm::mat4 mViewProjectionMatrix = pCamera->GetViewProjectionMatrix();
	GLint nLocationVPMatrix = glGetUniformLocation(pShader->GetID(), "ViewProjectionMatrix");
	glUniformMatrix4fv(nLocationVPMatrix, 1, GL_FALSE, &(mViewProjectionMatrix[0][0]));

	GLint nLocationObjectPosition = glGetUniformLocation(pShader->GetID(), "ObjectPosition");
	glUniform3f(nLocationObjectPosition, m_vParticlePosition.x, m_vParticlePosition.y, m_vParticlePosition.z);
	m_pModel->Draw(pShader);
	
	if (m_eDGL == EDGL::DoublePendulum)
	{
		glUniform3f(nLocationObjectPosition, m_vParticlePosition2.x, m_vParticlePosition2.y, m_vParticlePosition2.z);
		m_pModel->Draw(pShader);	
	}
	if (m_eDGL == EDGL::DoublePendulum || m_eDGL == EDGL::PendulumAngular || m_eDGL == EDGL::PendulumPosition || m_eDGL == EDGL::PendulumSpring)
	{
		glUniform3f(nLocationObjectPosition, 0, 0, 0);
		m_pModel->Draw(pShader);
	}
	pShader->UnBind();

	if (m_eDGL == EDGL::DoublePendulum)
	{
		pShaderLine->Bind();
		glBindVertexArray(m_nLineVAO);
		GLint nLocationObjectPosition = glGetUniformLocation(pShaderLine->GetID(), "ObjectPosition");
		glUniform3f(nLocationObjectPosition, m_vParticlePosition.x, m_vParticlePosition.y, m_vParticlePosition.z);
		GLint nLocationObjectPosition2 = glGetUniformLocation(pShaderLine->GetID(), "ObjectPosition2");
		glUniform3f(nLocationObjectPosition2, m_vParticlePosition2.x, m_vParticlePosition2.y, m_vParticlePosition2.z);

		
		glEnableVertexAttribArray(0);
		// pass information about how vertex array is composed
		glVertexAttribPointer(0, sizeof(glm::vec4), GL_FLOAT, GL_FALSE, sizeof(glm::vec4), NULL);

		glDrawArrays(GL_LINES, 0, 3);

		glBindVertexArray(0);
		pShaderLine->UnBind();
	}

}

std::vector<long double> Scale(std::vector<long double> vec, double scale)
{
	for (unsigned int i = 0; i < vec.size(); ++i)
	{
		vec[i] = vec[i] * scale;
	}
	return vec;
}

std::vector<long double> Add(std::vector<long double> vec1, std::vector<long double> vec2)
{
	for (unsigned int i = 0; i < vec1.size(); ++i)
	{
		vec1[i] = vec1[i] + vec2[i];
	}
	return vec1;
}

std::vector<long double> CPhysics::RungeKutta(std::vector<long double> state, double h)
{
	std::vector<long double> k1 = Scale(CalculateNewDerivation(state), h);
	std::vector<long double> k2 = Scale(CalculateNewDerivation(Add(state, Scale(k1, 0.5))), h);
	std::vector<long double> k3 = Scale(CalculateNewDerivation(Add(state, Scale(k2, 0.5))), h);
	std::vector<long double> k4 = Scale(CalculateNewDerivation(Add(state, k3)), h);

	std::vector<long double> k14 = Scale(Add(k1, k4), 1.0 / 6.0);
	std::vector<long double> k23 = Scale(Add(k2, k3), 1.0 / 3.0);

	return Add(Add(k14, k23), state);
}

float euler(float x, float dx, float h)
{
	return  x + h*dx;
}

void CPhysics::SimulateOneStep(float h)
{
	if(m_eDGL != EDGL::PendulumPosition)
	{

		if (true)
		{
			m_sParticle.vecVariables = RungeKutta(m_sParticle.vecVariables, h);
		}
		else
		{
			std::vector<long double> derivatives = CalculateNewDerivation(m_sParticle.vecVariables);
			for (unsigned int i = 0; i < m_sParticle.vecVariables.size(); ++i)
			{
				m_sParticle.vecVariables[i] = euler(m_sParticle.vecVariables[i], derivatives[i], h);
			}
		}		
	}
	else
	{
		glm::vec2 g = glm::vec2(0, 9.807); //gravitations konstante  "Pendulum"
		float l = 1;
		double m = 1;

		glm::vec2 vPositionOld = glm::vec2(m_sParticle.vecVariables[0], m_sParticle.vecVariables[1]);
		glm::vec2 vVelocityOld = glm::vec2(m_sParticle.vecVariables[2], m_sParticle.vecVariables[3]);
		glm::vec2 vAccelerationOld = glm::vec2(m_sParticle.vecVariables[4], m_sParticle.vecVariables[5]);

		glm::vec2 vCorrection = (glm::normalize(vPositionOld + (-g*h*h))*l-(vPositionOld + (-g*h*h)))/h;
		//glm::vec2 vAcceleration = (glm::normalize(vPositionOld +(-g*h*h)) - vPositionOld)/h;
		glm::vec2 vAcceleration = -g*h + vCorrection;

		glm::vec2 vVelocityNew = vVelocityOld + vAcceleration;

		glm::vec2 vPositionNew = vPositionOld + vVelocityNew*h;

		m_sParticle.vecVariables[0] = vPositionNew.x;
		m_sParticle.vecVariables[1] = vPositionNew.y;
		m_sParticle.vecVariables[2] = vVelocityNew.x;
		m_sParticle.vecVariables[3] = vVelocityNew.y;
		m_sParticle.vecVariables[4] = vAcceleration.x;
		m_sParticle.vecVariables[5] = vAcceleration.y;
		m_sParticle.vecVariables[6] = vCorrection.x;
		m_sParticle.vecVariables[7] = vCorrection.y;
	}
	this->SaveCurrentParticlePosition();
}

int CPhysics::GetCountTrackedTimeSteps()
{
	return m_vecRecordedData.size();
}















