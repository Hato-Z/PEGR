#version 120

////////////////////////////////////////////////////ADJUSTABLE VARIABLES/////////////////////////////////////////////////////////

#define POM 								//Comment to disable parallax occlusion mapping.
#define NORMAL_MAP_MAX_ANGLE 0.88f   		//The higher the value, the more extreme per-pixel normal mapping (bump mapping) will be.





/* Here, intervalMult might need to be tweaked per texture pack.  
   The first two numbers determine how many samples are taken per fragment.  They should always be the equal to eachother.
   The third number divided by one of the first two numbers is inversely proportional to the range of the height-map. */

const vec3 intervalMult = vec3(0.00058828125, 0.00058828125, 0.085); // Fine for 128x128 tile size

uniform sampler2D texture;
uniform sampler2D normals;
uniform sampler2D specular;

uniform float wetness;
uniform int worldTime;
uniform vec3 sunPosition;
uniform vec3 moonPosition;
varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 viewVector;

varying vec3 tangent;
varying vec3 normal;
varying vec3 binormal;

varying float translucent;
varying float distance;
varying float test;

const float MAX_OCCLUSION_DISTANCE = 100.0;

const int MAX_OCCLUSION_POINTS = 20;

const int GL_LINEAR = 9729;
const int GL_EXP = 2048;

uniform int fogMode;

float wetx = clamp(wetness, 0.0f, 1.0)/1.0;

const float bump_distance = 60.0;
const float pom_distance = 30.0;
const float fademult = 0.1;





void main() {	

	
	vec2 adjustedTexCoord = texcoord.st;
	float texinterval = 0.0625;
	float pomsample = 0.0;

#ifdef POM
	if (viewVector.z < 0.0 && distance < pom_distance && test < 0.2) {
		vec3 coord = vec3(texcoord.st, 1.0);

		if (texture2D(normals, coord.st).a < 1.0) {
			vec2 minCoord = vec2(texcoord.s - mod(texcoord.s, texinterval), texcoord.t - mod(texcoord.t, texinterval));
			vec2 maxCoord = vec2(minCoord.s + texinterval, minCoord.t + texinterval);
		
			vec3 interval = viewVector * intervalMult * 2.0;

			for (int loopCount = 0; pomsample < coord.z && loopCount < 14; ++loopCount) {
				coord += interval * clamp((1.0 - pomsample) * 10000.0f, 0.0, 1.0);
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
				pomsample = texture2D(normals, coord.st).a;
			}
		}
		else pomsample = 1.0;
		adjustedTexCoord = coord.st;
	}
	else pomsample = 1.0;
#endif


				  
	

	
	vec3 lightVector;
	if (worldTime < 12700 || worldTime > 23250) {
		lightVector = normalize(sunPosition);
	} else {
		lightVector = normalize(moonPosition);
	}
	vec3 indlmap = mix(min(lmcoord.t+0.1,1.0),1.0,lmcoord.s)*texture2D(texture,adjustedTexCoord).rgb*color.rgb;
	gl_FragData[0] = vec4(indlmap,texture2D(texture,adjustedTexCoord).a*color.a);
	
	
	
	
float pomdepth = 0.0;
	vec4 frag2 = vec4(vec3(normal) * 0.5 + 0.5, 1.0f);
	float dirtest;
	if (translucent > 0.9) dirtest = 0.4;
	else {
	dirtest = 1.0-0.8*step(dot(frag2.xyz*2.0-1.0,lightVector),-0.02);
	pomdepth = -pomsample*1.5-1.5;			//make the block closer as it is for self-shadowing fix
	}
	float pomdepthbias = (pomdepth) * (1.0 - gl_FragCoord.z) * (1.0 - gl_FragCoord.z);
	gl_FragDepth = gl_FragCoord.z+pomdepthbias;
	
	if (distance < bump_distance) {
	
			vec3 bump = texture2D(normals, adjustedTexCoord).rgb * 2.0 - 1.0;
			
			float bumpmult = clamp(bump_distance * fademult - distance * fademult, 0.0f, 1.0f) * NORMAL_MAP_MAX_ANGLE;
	
			bump = bump * vec3(bumpmult, bumpmult, bumpmult) + vec3(0.0f, 0.0f, 1.0f - bumpmult);
		mat3 tbnMatrix = mat3(tangent.x, binormal.x, normal.x,
								  tangent.y, binormal.y, normal.y,
						     	  tangent.z, binormal.z, normal.z);
			
			frag2 = vec4(bump * tbnMatrix * 0.5 + 0.5, 1.0);
			
	} 
	//normals
	gl_FragData[2] = frag2;	
	
	vec3 specularity = texture2D(specular,texcoord.xy).rgb;
			float g_spec = specularity.r + specularity.g* wetx;
			float g_irr = specularity.b;
			
			float totalspec = (g_spec + g_irr) * 0.333 * mix(min(lmcoord.t+0.1,1.0),1.0,lmcoord.s);
	//x = specularity / y = land(0.0/1.0)/shadow early exit(0.2)/water(0.05)/translucent(0.4) / z = torch lightmap
	
	gl_FragData[4] = vec4(totalspec, dirtest, lmcoord.s, 0.0);
	
}