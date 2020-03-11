#pragma once
#include <iostream>
#include <ctime>
#include "../ImGUI/imgui.h"
#include "../ImGUI/imgui_impl_glfw.h"
#include "../ImGUI/imgui_impl_opengl3.h"
#include <sstream>
#include "WindowGLFW.h"
#include "Camera.h"
#include "Shader.h"
#include "Physics.h"

CWindowGLFW* pWindow = nullptr;
CCamera* pCamera = nullptr;
CShader* pShader = nullptr;
CShader* pShaderLine = nullptr;
CPhysics* pPhysics = nullptr;

void LoadContent()
{
	pWindow = new CWindowGLFW(1200, 800);
	pCamera = new CCamera(glm::vec3(0, 0, -5), glm::vec3(0, 0, 1), glm::vec3(0, 1, 0), 110, static_cast<float>(pWindow->GetWindowSize().x) / static_cast<float>(pWindow->GetWindowSize().y), 0.01f, 1000.0f);
	pShader = new CShader();
	pShaderLine = new CShader();
	pPhysics = new CPhysics();
	pPhysics->SetInitialValues(0.01f);

	pShader->CreateShaderProgram("../shader/VS_showSphere.glsl", nullptr, nullptr, nullptr, "../shader/FS_showSphere.glsl");
	pShaderLine->CreateShaderProgram("../shader/VS_showLine.glsl", nullptr, nullptr, nullptr, "../shader/FS_showSphere.glsl");

	{ //GUI
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); (void)io;
		ImGui::StyleColorsDark();
		ImGui_ImplGlfw_InitForOpenGL(pWindow->GetWindowID(), true);
		ImGui_ImplOpenGL3_Init("#version 430");
	}
}

void InputManagement(double dElapsed_ms)
{
	float fTranslationSpeed = 0.1f;
	float fRotationSpeed = 0.01f;
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_W) == GLFW_PRESS)
	{
		pCamera->Translate(glm::vec3(0, 0, fTranslationSpeed));
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_S) == GLFW_PRESS)
	{
		pCamera->Translate(glm::vec3(0, 0, -fTranslationSpeed));
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_A) == GLFW_PRESS)
	{
		pCamera->Translate(glm::vec3(fTranslationSpeed, 0, 0));
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_D) == GLFW_PRESS)
	{
		pCamera->Translate(glm::vec3(-fTranslationSpeed, 0, 0));
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_LEFT) == GLFW_PRESS)
	{
		pCamera->Rotate(glm::vec3(0, 1, 0), -fRotationSpeed);
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_RIGHT) == GLFW_PRESS)
	{
		pCamera->Rotate(glm::vec3(0, 1, 0), fRotationSpeed);
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_UP) == GLFW_PRESS)
	{
		pCamera->Rotate(glm::vec3(1, 0, 0), -fRotationSpeed);
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_DOWN) == GLFW_PRESS)
	{
		pCamera->Rotate(glm::vec3(1, 0, 0), fRotationSpeed);
	}
	if (glfwGetKey(pWindow->GetWindowID(), GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{

	}
}

void GameLoop()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	bool bExitGame = false;
	clock_t cTimeStamp = clock();
	bool bEnableWireFrame = false;
	int nSimulationStepsPerFrame = 1;
	int nCountRecordedStates = 1000;
	double dSimulationTimeStepsSize = 0.01;
	bool bEnableSimulation = false;

	while (!bExitGame)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear the color and the depth buffer

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		bExitGame = pWindow->ManageInputs();

		clock_t cNewTimeStamp = clock();
		double dElapsed_ms = double(cNewTimeStamp - cTimeStamp) / (CLOCKS_PER_SEC * 1000.0);
		cTimeStamp = cNewTimeStamp;
		InputManagement(dElapsed_ms);


		if (bEnableWireFrame)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		{
			ImGui::Begin("Simulation");
		
			const char* items[] = { "PendulumAngular", "PendulumPosition", "PendulumSpring", "LorenzAttractor", "SinglePoint", "DoublePendulum"};
			static const char* current_item = NULL;

			if (ImGui::BeginCombo("DGL", current_item)) // The second parameter is the label previewed before opening the combo.
			{
				for (int n = 0; n < IM_ARRAYSIZE(items); n++)
				{
					bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
					if (ImGui::Selectable(items[n], is_selected))
						current_item = items[n];
					if (is_selected)
						ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
				}
				ImGui::EndCombo();
			}
			pPhysics->SetDGL(current_item);

			if (ImGui::Button("Start simulation"))
			{
				pPhysics->SetInitialValues(dSimulationTimeStepsSize);
				bEnableSimulation = true;
			}
			if (ImGui::Button("Stop simulation"))
			{
				bEnableSimulation = false;
			}

			ImGui::InputInt("Count Simulation-\nSteps per Frame", &nSimulationStepsPerFrame);
			if (nSimulationStepsPerFrame < 1)
			{
				nSimulationStepsPerFrame = 1;
			}
			ImGui::InputDouble("Simulation Time-\nStep Size", &dSimulationTimeStepsSize);
		
			double dL1 = pPhysics->GetL1();
			double dL2 = pPhysics->GetL2();
			double dM1 = pPhysics->GetM1();
			double dM2 = pPhysics->GetM2();
			double dK = pPhysics->GetK();

			if(current_item == "PendulumAngular" || current_item == "PendulumPosition" || current_item == "PendulumSpring"|| current_item == "DoublePendulum")
				ImGui::InputDouble("Length1", &dL1);
			if (current_item == "DoublePendulum")
				ImGui::InputDouble("Length2", &dL2);
			if (current_item == "DoublePendulum")
				ImGui::InputDouble("Mass1", &dM1);
			if (current_item == "DoublePendulum")
				ImGui::InputDouble("Mass2", &dM2);
			if (current_item == "PendulumSpring")
				ImGui::InputDouble("Spring constant", &dK);

			pPhysics->SetL1(dL1);
			pPhysics->SetL2(dL2);
			pPhysics->SetM1(dM1);
			pPhysics->SetM2(dM2);
			pPhysics->SetK(dK);

			ImGui::InputInt("#states to record", &nCountRecordedStates);

			pPhysics->SetCountRecordedStates(nCountRecordedStates);
			ImGui::End();

			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		}

		//RunPhysics

		if (bEnableSimulation == true && nCountRecordedStates <= pPhysics->GetCountTrackedTimeSteps() && nCountRecordedStates > 0)
		{
			bEnableSimulation = false;

			pPhysics->SaveTrackedValues();
		}
		if (bEnableSimulation)
		{
			pPhysics->RecordState();
			for (unsigned int i = 0; i < nSimulationStepsPerFrame; ++i)
			{
				pPhysics->SimulateOneStep(dSimulationTimeStepsSize);
			}
		}
		//else
		{
			pPhysics->DrawParticle(pShader, pShaderLine, pCamera);
			pWindow->SwapBuffers();
		}
	}

}

void DeleteContent()
{
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();
	delete pWindow;
	delete pCamera;
	delete pShader;
	delete pShaderLine;
	delete pPhysics;
}

void main()
{
	LoadContent();
	GameLoop();
	DeleteContent();
}
