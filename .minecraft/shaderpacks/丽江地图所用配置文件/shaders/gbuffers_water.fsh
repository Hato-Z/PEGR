#version 120

uniform sampler2D texture;
uniform sampler2D lightmap;

uniform mat4 gbufferProjectionInverse;
uniform mat4 gbufferPreviousProjection;

uniform mat4 gbufferModelViewInverse;
uniform mat4 gbufferPreviousModelView;

uniform vec3 cameraPosition;
uniform vec3 previousCameraPosition;

uniform float rainStrength;

varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying float iswater;
varying float isice;

//varying vec4 bloommask;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;

void main() {

	vec4 tex = texture2D(texture, texcoord.st);
	
	float zero = 1.0f;
	float transx = 0.0f;
	float transy = 0.0f;
	//float iswater = 0.0f;
	
	float texblock = 0.0625f;

	
	if (iswater > 0.999f) {
		tex = vec4(0.9f, 0.9f, 0.9f, 0.15f);
	}
	

/*
	if (texcoord.s >= 0.8125f && texcoord.t >= 0.75f && texcoord.t <= 0.875f) {
		tex = vec4(0.5f, 0.9f, 0.9f, 0.2f);
		iswater = 1.0f;
	}	else {
		iswater = 0.0f;
	}
*/


	gl_FragData[0] = tex * color;
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z), 1.0);
	gl_FragData[2] = vec4(vec3(normal.x, normal.y, normal.z) * 0.5 + 0.5, 1.0f);
	gl_FragData[4] = vec4(0.0, iswater, 1.0, 1.0);
	//store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	gl_FragData[5] = vec4(texture2D(lightmap, lmcoord.st).rgb, 0.0f);
	gl_FragData[6] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

	float fogsat = 1.3;
	
	vec3 fogcolor = gl_Fog.color.rgb;
	
	fogcolor.r = mix(fogcolor.r * 0.8, fogcolor.r * 0.9, rainx);
	fogcolor.g = mix(fogcolor.g * 0.8, fogcolor.g * 0.9, rainx);
	fogcolor.b = mix(fogcolor.b * 0.8, fogcolor.b * 0.9, rainx);
	
	fogsat = mix(1.3, 0.8, rainx);

	

	
	fogcolor.r = (fogcolor.r * fogsat) - (((fogcolor.g + fogcolor.b) / 2.0) * (fogsat - 1.0));
	fogcolor.g = (fogcolor.g * fogsat) - (((fogcolor.r + fogcolor.b) / 2.0) * (fogsat - 1.0));
	fogcolor.b = (fogcolor.b * fogsat) - (((fogcolor.r + fogcolor.g) / 2.0) * (fogsat - 1.0));
	
	
	//gl_FragData[1] = vec4(0.0);
		
		
	if (fogMode == GL_EXP) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, fogcolor.rgb, 1.0 - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
	} else if (fogMode == GL_LINEAR) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, fogcolor.rgb, clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
	}
}