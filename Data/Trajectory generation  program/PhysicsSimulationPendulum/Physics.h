#pragma once
#include <vector>
#include <string>
#include "Shader.h"
#include "Camera.h"
#include "ModelObject.h"
#include <glad/glad.h>

class CPhysics
{
public:
	CPhysics();
	~CPhysics();

	struct SState
	{
		std::vector<long double> vecVariables;
	};

	void SimulateOneStep(float h);
	void SetInitialValues(float h);
	void DrawParticle(CShader* pShader, CShader* pShaderLine, CCamera* pCamera);
	void SaveTrackedValues();
	void SetCountRecordedStates(int nCount);
	void RecordState();
	void SaveCurrentParticlePosition();
	int GetCountTrackedTimeSteps();
	void SetDGL(const char* pDGL);
	void SetL1(double l1);
	void SetL2(double l2);
	void SetM1(double m1);
	void SetM2(double m2);
	void SetK(double k);

	double GetL1();
	double GetL2();
	double GetM1();
	double GetM2();
	double GetK();
	enum EDGL
	{
		PendulumAngular,
		PendulumPosition,
		PendulumSpring,
		LorenzAttractor,
		SinglePoint,
		EasyExample,
		DoublePendulum
	};

private:
	std::vector<long double> RungeKutta(std::vector<long double> state, double h);
	std::vector<long double> CalculateNewDerivation(std::vector<long double> particle);

	SState m_sParticle;
	int m_nCountRecordedStates = -1;
	std::string m_strFileName;
	glm::vec3 m_vParticlePosition;
	glm::vec3 m_vParticlePosition2;
	std::vector<long double> m_vAcceleration;
	std::vector<SState> m_vecRecordedData;
	CModel* m_pModel;
	GLuint m_nLineVBO;
	GLuint m_nLineVAO;
	double m_dL1 = 1;
	double m_dL2 = 1;
	double m_dM1 = 1;
	double m_dM2 = 1;
	double m_dK = 1000;

	EDGL m_eDGL;
};

