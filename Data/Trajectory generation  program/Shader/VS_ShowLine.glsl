#version 430 core

layout(location = 0) in vec3 PositionVertex;
uniform mat4 ViewProjectionMatrix;
uniform vec3 ObjectPosition;
uniform vec3 ObjectPosition2;

void main()
{	
	int vID = gl_VertexID;
	if(vID == 1)
	{
		gl_Position = ViewProjectionMatrix * vec4(ObjectPosition, 1.0f);
	}
	else
	{
		gl_Position = ViewProjectionMatrix * vec4(ObjectPosition2, 1.0f);
	}
}