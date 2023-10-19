#version 120

varying vec4 texcoord;

uniform int worldTime;

uniform vec3 sunPosition;
uniform vec3 moonPosition;

varying vec3 lightVector;

void main() {
	gl_Position = ftransform();
	
	if (worldTime < 12700 || worldTime > 23250) {
		lightVector = normalize(sunPosition);
	} else {
		lightVector = normalize(moonPosition);
	}
	
	texcoord = gl_MultiTexCoord0;
}
