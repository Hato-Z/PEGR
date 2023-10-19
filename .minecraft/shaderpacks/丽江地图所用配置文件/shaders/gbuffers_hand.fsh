#version 120


////////////////////////////////////////////////////ADJUSTABLE VARIABLES/////////////////////////////////////////////////////////

#define POM 								//Comment to disable parallax occlusion mapping.
#define NORMAL_MAP_MAX_ANGLE 0.88f   		//The higher the value, the more extreme per-pixel normal mapping (bump mapping) will be.





/* Here, intervalMult might need to be tweaked per texture pack.  
   The first two numbers determine how many samples are taken per fragment.  They should always be the equal to eachother.
   The third number divided by one of the first two numbers is inversely proportional to the range of the height-map. */

//const vec3 intervalMult = vec3(0.0039, 0.0039, 4.5); // Fine for 16x16 tile size
//const vec3 intervalMult = vec3(0.0019, 0.0019, 0.5); // Fine for 32x32 tile size
//const vec3 intervalMult = vec3(0.00048828125, 0.00048828125, 0.2); // Fine for 128x128 tile size
const vec3 intervalMult = vec3(0.00058828125, 0.00058828125, 0.085); // Fine for 128x128 tile size

uniform sampler2D texture;
uniform sampler2D lightmap;
uniform sampler2D normals;
uniform sampler2D specular;
uniform float rainStrength;


varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 viewVector;
varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;

varying float translucent;
varying float distance;

const float MAX_OCCLUSION_DISTANCE = 100.0;

const int MAX_OCCLUSION_POINTS = 20;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;

const float bump_distance = 80.0f;
const float fademult = 0.1f;





void main() {	



		
	
	vec2 adjustedTexCoord = texcoord.st;
	float texinterval = 0.0625f;

#ifdef POM
	if (viewVector.z < 0.0) {
		vec3 coord = vec3(texcoord.st, 1.0);

		if (texture2D(normals, coord.st).a < 1.0) {
			vec2 minCoord = vec2(texcoord.s - mod(texcoord.s, texinterval), texcoord.t - mod(texcoord.t, texinterval));
			vec2 maxCoord = vec2(minCoord.s + texinterval, minCoord.t + texinterval);
		
			vec3 interval = viewVector * intervalMult * 2.0f;

			for (int loopCount = 0; texture2D(normals, coord.st).a < coord.z && loopCount < 14; ++loopCount) {
				coord += interval * clamp((1.0f - texture2D(normals, coord.st).a) * 10000.0f, 0.0f, 1.0f);
				if (coord.s < minCoord.s) {
					coord.s += texinterval;
				} else if (coord.s >= maxCoord.s) {
					coord.s -= texinterval;
				}
				if (coord.t < minCoord.t) {
					coord.t += texinterval;
				} else if (coord.t >= maxCoord.t) {
					coord.t -= texinterval;
				}
			}
		}

		adjustedTexCoord = coord.st;
	}
#endif

	float alpha = texture2D(texture, adjustedTexCoord).a;
	float alphacheck = 0.0f;
	
		if (alpha > 0.9) {
			alphacheck = 1.0f;
		}
				  
	float pomdepth = texture2D(normals, adjustedTexCoord).a;

	float pomdepthbias = (1.0f - pomdepth) * (1.0f - gl_FragCoord.z) * (1.0f - gl_FragCoord.z);

	gl_FragData[0] = texture2D(texture, adjustedTexCoord) * (color);
	gl_FragData[1] = vec4(vec3(gl_FragCoord.z) + vec3(pomdepthbias) * 2.0f, alphacheck);
	
	
	vec4 frag2;
	
	if (distance < bump_distance) {
	
			vec3 bump = texture2D(normals, adjustedTexCoord).rgb * 2.0f - 1.0f;
			
			float bumpmult = clamp(bump_distance * fademult - distance * fademult, 0.0f, 1.0f) * NORMAL_MAP_MAX_ANGLE;
	
			bump = bump * vec3(bumpmult, bumpmult, bumpmult) + vec3(0.0f, 0.0f, 1.0f - bumpmult);
		
			mat3 tbnMatrix = mat3(tangent.x, binormal.x, normal.x,
								  tangent.y, binormal.y, normal.y,
						     	  tangent.z, binormal.z, normal.z);
			
			frag2 = vec4(bump * tbnMatrix * 0.5 + 0.5, 1.0);
			
	} else {
	
			frag2 = vec4((normal) * 0.5f + 0.5f, 1.0f);		
			
	}
	
	gl_FragData[2] = frag2;	
	
	gl_FragData[4] = vec4(translucent, 0.0f, 1.0f, 0.0f);
	
	//store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
	
	vec4 lightmap = vec4(texture2D(lightmap, lmcoord.st).rgb, alphacheck);
		
		if (fogMode == GL_EXP) {
			lightmap = mix(lightmap, vec4(0.0f, 0.0f, 1.0f, 1.0f),  1.0f - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
		} else if (fogMode == GL_LINEAR) {
			lightmap = mix(lightmap, vec4(0.0f, 0.0f, 1.0f, 1.0f),  clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
		}
	
	
	gl_FragData[5] = lightmap;
	gl_FragData[6] = texture2D(specular, adjustedTexCoord.st);
	
	
	float fogsat = 1.3;
	
	vec3 fogcolor = gl_Fog.color.rgb;
	
	fogcolor.r = mix(fogcolor.r * 0.5f, fogcolor.r * 0.9f, rainx);
	fogcolor.g = mix(fogcolor.g * 0.5f, fogcolor.g * 0.9f, rainx);
	fogcolor.b = mix(fogcolor.b * 0.5f, fogcolor.b * 0.9f, rainx);
	
	fogsat = mix(1.3f, 0.8f, rainx);

	

	
	fogcolor.r = (fogcolor.r * fogsat) - (((fogcolor.g + fogcolor.b) / 2.0) * (fogsat - 1.0));
	fogcolor.g = (fogcolor.g * fogsat) - (((fogcolor.r + fogcolor.b) / 2.0) * (fogsat - 1.0));
	fogcolor.b = (fogcolor.b * fogsat) - (((fogcolor.r + fogcolor.g) / 2.0) * (fogsat - 1.0));
	
	

	if (fogMode == GL_EXP) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (fogcolor.rgb), 1.0 - clamp(exp(-gl_Fog.density * gl_FogFragCoord), 0.0, 1.0));
	} else if (fogMode == GL_LINEAR) {
		gl_FragData[0].rgb = mix(gl_FragData[0].rgb, (fogcolor.rgb), clamp((gl_FogFragCoord - gl_Fog.start) * gl_Fog.scale, 0.0, 1.0));
	}
}