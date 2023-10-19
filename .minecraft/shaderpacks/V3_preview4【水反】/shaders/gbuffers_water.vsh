#version 120

#define WAVING_WATER

uniform int worldTime;

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 normal;

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
	vec4 positionLeft  = position + vec4(-1.0f, 0.0f,  0.0f, 0.0f);
	vec4 positionRight = position + vec4( 1.0f, 0.0f,  0.0f, 0.0f);
	vec4 positionUp    = position + vec4( 0.0f, 0.0f,  1.0f, 0.0f);
	vec4 positionDown  = position + vec4( 0.0f, 0.0f, -1.0f, 0.0f);
	float displaceXC = 0.0;
	float displaceYC = 0.0;
	float displaceXL = 0.0;
	float displaceYL = 0.0;
	float displaceXU = 0.0;
	float displaceYU = 0.0;
	
	float octaveScale   = 0.05f;
	float waveAmplitude = 1.18f;
	float waveSpeed 	= 0.81f;
	
	vec2 displaceDelta = vec2(0.0f);
	
	
		if (mc_Entity.x == 8 || mc_Entity.x == 9) {
			iswater = 1.0f;
			
			#ifdef WAVING_WATER
			//Center displacement
			for(int i = 1; i < 5; ++i){
			
				float octave = pow(2.0f, i) * octaveScale;
				float speed = i * waveSpeed;
				
				float magnitude = (sin((position.y * octave + position.x * octave + worldTime * octave * 3.14159265358979323846264 / ((28.0) * speed))) * 0.15 + 0.15) * waveAmplitude / pow(2.0f, 6 - i);
				float d0 = sin(position.y * octave * 3.0 + position.x * octave * 0.3 + worldTime * 3.14159265358979323846264 / (112.0 * speed)) * 3.0 - 1.5;
				float d1 = sin(position.y * octave * 0.7 - position.x * octave * 10.0 + worldTime * 3.14159265358979323846264 / (142.0 * speed)) * 3.0 - 1.5;
				float d2 = sin(worldTime * 3.14159265358979323846264 / (132.0 * speed)) * 3.0 - 1.5;
				float d3 = sin(worldTime * 3.14159265358979323846264 / (122.0 * speed)) * 3.0 - 1.5;
				
				displaceXC += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave + d2) + (position.x * octave + d3)) * (magnitude);
				displaceYC += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave * -0.5 + d1) + (position.x * octave * -2.0 + d0)) * (magnitude);
				
				position.y += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave + d2) + (position.x * octave + d3)) * (magnitude);
				position.y -= sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (position.z * octave * -0.5 + d1) + (position.x * octave * -2.0 + d0)) * (magnitude);
				
			}
			
			//Left displacement
			for(int i = 1; i < 5; ++i){
			
				float octave = pow(2.0f, i) * octaveScale;
				float speed = i * waveSpeed;
				
				float magnitude = (sin((positionLeft.y * octave + positionLeft.x * octave + worldTime * octave * 3.14159265358979323846264 / ((28.0) * speed))) * 0.15 + 0.15) * waveAmplitude / pow(2.0f, 6 - i);
				float d0 = sin(positionLeft.y * octave * 3.0 + positionLeft.x * octave * 0.3 + worldTime * 3.14159265358979323846264 / (112.0 * speed)) * 3.0 - 1.5;
				float d1 = sin(positionLeft.y * octave * 0.7 - positionLeft.x * octave * 10.0 + worldTime * 3.14159265358979323846264 / (142.0 * speed)) * 3.0 - 1.5;
				float d2 = sin(worldTime * 3.14159265358979323846264 / (132.0 * speed)) * 3.0 - 1.5;
				float d3 = sin(worldTime * 3.14159265358979323846264 / (122.0 * speed)) * 3.0 - 1.5;
				
				displaceXL += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (positionLeft.z * octave + d2) + (positionLeft.x * octave + d3)) * (magnitude);
				displaceYL += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (positionLeft.z * octave * -0.5 + d1) + (positionLeft.x * octave * -2.0 + d0)) * (magnitude);
				
			}
			
			//Up displacement
			for(int i = 1; i < 5; ++i){
			
				float octave = pow(2.0f, i) * octaveScale;
				float speed = i * waveSpeed;
				
				float magnitude = (sin((positionUp.y * octave + positionUp.x * octave + worldTime * octave * 3.14159265358979323846264 / ((28.0) * speed))) * 0.15 + 0.15) * waveAmplitude / pow(2.0f, 6 - i);
				float d0 = sin(positionUp.y * octave * 3.0 + positionUp.x * octave * 0.3 + worldTime * 3.14159265358979323846264 / (112.0 * speed)) * 3.0 - 1.5;
				float d1 = sin(positionUp.y * octave * 0.7 - positionUp.x * octave * 10.0 + worldTime * 3.14159265358979323846264 / (142.0 * speed)) * 3.0 - 1.5;
				float d2 = sin(worldTime * 3.14159265358979323846264 / (132.0 * speed)) * 3.0 - 1.5;
				float d3 = sin(worldTime * 3.14159265358979323846264 / (122.0 * speed)) * 3.0 - 1.5;
				
				displaceXU += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (positionUp.z * octave + d2) + (positionUp.x * octave + d3)) * (magnitude);
				displaceYU += sin((worldTime * 3.14159265358979323846264 / (11.0 * speed)) + (positionUp.z * octave * -0.5 + d1) + (positionUp.x * octave * -2.0 + d0)) * (magnitude);
				
			}
			
			displaceDelta = vec2(displaceXC - displaceXL, displaceYU - displaceYC);
			
			
			#endif
		}
	
		


	gl_Position = gl_ProjectionMatrix * (gl_ModelViewMatrix * position);

	
	color = gl_Color;
	
	texcoord = gl_TextureMatrix[0] * gl_MultiTexCoord0;

	lmcoord = gl_TextureMatrix[1] * gl_MultiTexCoord1;
	


	gl_FogFragCoord = gl_Position.z;
	
	
	normal = normalize(gl_NormalMatrix * gl_Normal );

	
}