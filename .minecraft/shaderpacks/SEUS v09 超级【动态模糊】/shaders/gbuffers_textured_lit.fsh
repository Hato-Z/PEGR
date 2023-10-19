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

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;

void main() {

	gl_FragData[0] = texture2D(texture, texcoord.st) * texture2D(lightmap, lmcoord.st) * (color*0.8);
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z), 1.0);
	gl_FragData[4] = vec4(0.0, 0.0, 1.0, 1.0);
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