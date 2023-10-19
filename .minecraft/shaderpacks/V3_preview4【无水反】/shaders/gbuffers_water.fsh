#version 120
uniform sampler2D texture;

uniform float rainStrength;

varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;


const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;

void main() {

	vec4 tex = texture2D(texture, texcoord.st);





	gl_FragData[0] = texture2D(texture,texcoord.xy)*color;
	gl_FragDepth = gl_FragCoord.z;
	gl_FragData[2] = vec4(normal * 0.5 + 0.5, 1.0f);
	//x = specularity / y = land(0.0/1.0)/shadow early exit(0.2)/water(0.05) / z = torch lightmap
	gl_FragData[4] = vec4(0.0, 1.0, lmcoord.s, 1.0);
}