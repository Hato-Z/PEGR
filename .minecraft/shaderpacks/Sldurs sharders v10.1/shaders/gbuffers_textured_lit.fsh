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
	float iswater = 0.0f;
	
	float texblock = 0.0625f;

	
	if (tex.r > 0.99 && tex.b > 0.99 && tex.g < 0.1f) {
		tex = vec4(0.5f, 0.9f, 0.9f, 0.2f);
		iswater = 1.0f;
	}	else {
		iswater = 0.0f;
	}


/*
	if (texcoord.s >= 0.8125f && texcoord.t >= 0.75f && texcoord.t <= 0.875f) {
		tex = vec4(0.5f, 0.9f, 0.9f, 0.2f);
		iswater = 1.0f;
	}	else {
		iswater = 0.0f;
	}
*/
	//store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	
	vec3 lightmaptorch = texture2D(lightmap, vec2(lmcoord.s, 0.00f)).rgb;
	vec3 lightmapsky   = texture2D(lightmap, vec2(0.0f, lmcoord.t)).rgb;
	
	//vec4 lightmap = texture2D(lightmap, lmcoord.st);
	vec4 lightmap = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	
		//Separate lightmap types
		lightmap.r = max(lightmaptorch.r, max(lightmaptorch.g, lightmaptorch.b)) * 2.0f;
		lightmap.b = lightmapsky.b;
		

	gl_FragData[0] = tex * color;
	gl_FragData[1] = vec4(gl_FragCoord.z, lightmap.r, lightmap.b, 1.0);
	
	
	//Material ID's
	
		//Scale material masks
		float land  	 = 1.0f/255.0f;

		
		//Combine material masks
		float mats_1 = land;	
	
	
		

	gl_FragData[2] = vec4(normal.rgb * 0.5f + 0.5f, 0.0f);

	
	
	gl_FragData[3] = vec4(mats_1, 0.0f, 0.0f, 1.0f);
	//gl_FragData[3] = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	//gl_FragData[6] = vec4(0.0f, 0.0f, 0.0f, 1.0f);


}