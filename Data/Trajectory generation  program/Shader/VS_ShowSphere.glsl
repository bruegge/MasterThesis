#version 430 core

layout(location = 0) in vec3 PositionVertex;
uniform mat4 ViewProjectionMatrix;
uniform vec3 ObjectPosition;

void main()
{	
	gl_Position = ViewProjectionMatrix * vec4(PositionVertex + ObjectPosition, 1.0f);
}