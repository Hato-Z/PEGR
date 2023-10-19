#version 120

uniform int worldTime;

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;

attribute vec4 mc_Entity;

varying float iswater;
varying float isice;



//varying vec4 bloommask;

//attribute vec4 mc_Entity;

void main() {

	iswater = 0.0f;
	isice = 0.0f;

	//bloommask = vec4(0.0);
	
	//if (mc_Entity.x == 8.0 || mc_Entity.x == 9.0) {
	//	bloommask.x = 1.0f;
	//}
	
	

	
	if (mc_Entity.x == 79) {
		isice = 1.0f;
	}
	
	
	vec4 position = gl_Vertex;
	
	if (mc_Entity.x == 8 || mc_Entity.x == 9) {
		iswater = 1.0f;
	
		for(int i = 1; i < 4; ++i){
		
			float octave = i * 0.5;
			float speed = (octave) * 2.0;
			
			float magnitude = (sin((position.y * octave + position.x * octave + worldTime * octave * 3.14159265358979323846264 / ((28.0) * speed))) * 0.15 + 0.15) * 0.28;
			float d0 = sin(position.y * octave * 3.0 + position.x * octave * 0.3 + worldTime * 3.14159265358979323846264 / (112.0 * speed)) * 3.0 - 1.5;
			float d1 = sin(position.y * octave * 0.7 - position.x * octave * 10.0 + worldTime * 3.14159265358979323846264 / (142.0 * speed)) * 3.0 - 1.5;
			float d2 = sin(worldTime * 3.14159265358979323846264 / (132.0 * speed)) * 3.0 - 1.5;
			float d3 = sin(worldTime * 3.14159265358979323846264 / (122.0 * speed)) * 3.0 - 1.5;
			position.y += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave + d2) + (position.x * octave + d3)) * (magnitude/2.0);
			position.y -= sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave * 0.5 + d1) + (position.x * octave * 2.0 + d0)) * (magnitude/2.0);
			
		}
	}
	vec4 locposition = gl_ModelViewMatrix * gl_Vertex;
	


	gl_Position = gl_ProjectionMatrix * (gl_ModelViewMatrix * position);

	
	color = gl_Color;
	
	texcoord = gl_TextureMatrix[0] * gl_MultiTexCoord0;

	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;
	


	gl_FogFragCoord = gl_Position.z;
	
	
	normal = normalize(gl_NormalMatrix * gl_Normal);

	if (gl_Normal.x > 0.5) {
		//  1.0,  0.0,  0.0
		tangent  = normalize(gl_NormalMatrix * vec3( 0.0,  0.0, -1.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.x < -0.5) {
		// -1.0,  0.0,  0.0
		tangent  = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.y > 0.5) {
		//  0.0,  1.0,  0.0
		tangent  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
	} else if (gl_Normal.y < -0.5) {
		//  0.0, -1.0,  0.0
		tangent  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0,  0.0,  1.0));
	} else if (gl_Normal.z > 0.5) {
		//  0.0,  0.0,  1.0
		tangent  = normalize(gl_NormalMatrix * vec3( 1.0,  0.0,  0.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	} else if (gl_Normal.z < -0.5) {
		//  0.0,  0.0, -1.0
		tangent  = normalize(gl_NormalMatrix * vec3(-1.0,  0.0,  0.0));
		binormal = normalize(gl_NormalMatrix * vec3( 0.0, -1.0,  0.0));
	}
	
}