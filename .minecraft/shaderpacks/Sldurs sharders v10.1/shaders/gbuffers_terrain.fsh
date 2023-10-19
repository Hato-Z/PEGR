#version 120


////////////////////////////////////////////////////ADJUSTABLE VARIABLES/////////////////////////////////////////////////////////

//#define POM 								//Comment to disable parallax occlusion mapping.
#define POM_HQ								//Higher quality parallax occlusion mapping.
#define NORMAL_MAP_MAX_ANGLE 0.25f   		//The higher the value, the more extreme per-pixel normal mapping (bump mapping) will be.





/* Here, intervalMult might need to be tweaked per texture pack.  
   The first two numbers determine how many samples are taken per fragment.  They should always be the equal to eachother.
   The third number divided by one of the first two numbers is inversely proportional to the range of the height-map. */

//const vec3 intervalMult = vec3(0.0039, 0.0039, 4.5); // Fine for 16x16 tile size
//const vec3 intervalMult = vec3(0.0019, 0.0019, 0.5); // Fine for 32x32 tile size
//const vec3 intervalMult = vec3(0.00048828125, 0.00048828125, 0.2); // Fine for 128x128 tile size

#ifdef POM_HQ
	const vec3 intervalMult = vec3(0.00015, 0.00015, 0.025); // Fine for 128x128 tile size
#else
	const vec3 intervalMult = vec3(0.00058828125, 0.00058828125, 0.085); // Fine for 128x128 tile size
#endif



uniform sampler2D texture;
uniform sampler2D lightmap;
uniform sampler2D normals;
uniform sampler2D specular;
//uniform float rainStrength;
uniform float wetness;


varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 viewVector;
varying vec3 normal;
varying vec3 tangent;
varying vec3 binormal;
varying vec3 globalNormal;

varying float leaves;
varying float grass;
varying float translucent;

varying float distance;

const float MAX_OCCLUSION_DISTANCE = 100.0;

const int MAX_OCCLUSION_POINTS = 20;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

const float bump_distance = 80.0f;
const float fademult = 0.1f;





void main() {	

if (texture2D(texture, texcoord.st).a == 0.0f){
	discard;
}


	//store lightmap in auxilliary texture. r = torch light. g = lightning. b = sky light.
		vec3 lightmaptorch = texture2D(lightmap, vec2(lmcoord.s, 0.00f)).rgb;
		vec3 lightmapsky   = texture2D(lightmap, vec2(0.0f, lmcoord.t)).rgb;
		
		//vec4 lightmap = texture2D(lightmap, lmcoord.st);
		vec4 lightmap = vec4(0.0f, 0.0f, 0.0f, 1.0f);
		
		//Separate lightmap types
		lightmap.r = max(lightmaptorch.r, max(lightmaptorch.g, lightmaptorch.b)) * 2.0f;
		lightmap.b = lightmapsky.b;

	
	vec2 adjustedTexCoord = texcoord.st;
	float texinterval = 0.0625f;
	
	vec3 coord = vec3(texcoord.st, 1.0);
	
#ifdef POM
	if (viewVector.z < 0.0 && distance < bump_distance && texture2D(normals, coord.st).a < 1.0f) {
				
		vec2 minCoord = vec2(texcoord.s - mod(texcoord.s, texinterval), texcoord.t - mod(texcoord.t, texinterval));
		vec2 maxCoord = vec2(minCoord.s + texinterval, minCoord.t + texinterval);
	
		vec3 interval = viewVector * intervalMult * 2.0f;
		
		#ifdef POM_HQ
		for (int loopCount = 0; texture2D(normals, coord.st).a < coord.z && loopCount < 40; ++loopCount) {
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
		
		#else
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
		#endif
	
		adjustedTexCoord = coord.st;
	}
#endif

				  
	float pomdepth = texture2D(normals, adjustedTexCoord).a;

	float pomdepthbias = (1.0f - pomdepth) * (1.0f - gl_FragCoord.z) * (1.0f - gl_FragCoord.z);

	//Diffuse
	gl_FragData[0] = texture2D(texture, adjustedTexCoord) * min((color * 1.5f), vec4(1.0f));
	
	//Depth
	gl_FragData[1] = vec4(gl_FragCoord.z + pomdepthbias * 1.5f, lightmap.r, lightmap.b, 1.0f);
	
	vec4 spec = texture2D(specular, adjustedTexCoord.st);
	
	
	vec4 frag2;
	vec4 frag3;
	
	if (distance < bump_distance) {
	
			vec3 bump = texture2D(normals, adjustedTexCoord).rgb * 2.0f - 1.0f;
			vec3 bump2 = texture2D(normals, adjustedTexCoord).rgb * 2.0f - 1.0f;
				 bump2.g = 1.0f - bump2.g;
			
			float bumpmult = clamp(bump_distance * fademult - distance * fademult, 0.0f, 1.0f) * NORMAL_MAP_MAX_ANGLE;
	              bumpmult *= 1.0f - (spec.g * 0.7f * wetness);
				  
			bump = bump * vec3(bumpmult, bumpmult, bumpmult) + vec3(0.0f, 0.0f, 1.0f - bumpmult);
			
			
			bump2 = bump2 * vec3(bumpmult, bumpmult, bumpmult) + vec3(0.0f, 0.0f, 1.0f - bumpmult);
		
			mat3 tbnMatrix = mat3(tangent.x, binormal.x, normal.x,
								  tangent.y, binormal.y, normal.y,
						     	  tangent.z, binormal.z, normal.z);
			
			frag2 = vec4(bump * tbnMatrix * 0.5 + 0.5, 1.0);
			
			frag3 = vec4((globalNormal + (bump2 * NORMAL_MAP_MAX_ANGLE)) * 0.5f + 0.5f, 1.0f);
			
	} else {
	
			frag2 = vec4((normal) * 0.5f + 0.5f, 1.0f);		
			frag3 = vec4(globalNormal * 0.5f + 0.5f, 1.0f);
			
	}
	
	
	
	

	
	//Material ID's
	
		//Scale material masks
		float land  	 = 1;

		
		//Combine material masks
		float mats_1 = land;
			  mats_1 = max(mats_1, grass);
			  mats_1 = max(mats_1, leaves);
			  mats_1 = max(mats_1, translucent);

	
	gl_FragData[2] = frag2;
		
	//matIDs, specularity, iswater
	gl_FragData[3] = vec4(mats_1/255.0f, spec.r + spec.b + spec.g * wetness, 0.0f, 1.0f);	
	
	//specularity, blank, iswater
	//gl_FragData[3] = vec4(spec.r + spec.b + spec.g * wetness, 0.0f, 0.0f, 1.0f);
	
	//gaux3
	//gl_FragData[6] = vec4(0.0f, 0.0f, 0.0f, 1.0f);

	

}