#version 120






/* Here, intervalMult might need to be tweaked per texture pack.  
   The first two numbers determine how many samples are taken per fragment.  They should always be the equal to eachother.
   The third number divided by one of the first two numbers is inversely proportional to the range of the height-map. */


uniform sampler2D texture;
uniform float rainStrength;

varying vec4 color;
varying vec4 texcoord;
varying vec4 lmcoord;

varying vec3 normal;

varying float iswater;


const float MAX_OCCLUSION_DISTANCE = 100.0;

const int MAX_OCCLUSION_POINTS = 20;

uniform int worldTime;

float rainx = clamp(rainStrength, 0.0f, 1.0f)/1.0f;

const float bump_distance = 80.0f;
const float fademult = 0.1f;





void main() {	
	

    vec4 tex = texture2D(texture, texcoord.xy);
    /*
    if (iswater > 0.9) {
		tex = vec4(	0.9,0.9,0.9, 0.15);
	}
    */

	
	
	vec3 indlmap = mix(lmcoord.t,1.0,lmcoord.s)*texture2D(texture,texcoord.xy).rgb*color.rgb;
	gl_FragData[0] = vec4(indlmap,texture2D(texture,texcoord.xy).a*color.a*0.85);
	gl_FragDepth = gl_FragCoord.z;
	
	vec4 frag2;
	

			frag2 = vec4((normal) * 0.5f + 0.5f, 1.0f);			
	
	
	gl_FragData[2] = frag2;	
	//x = specularity / y = land(0.0/1.0)/shadow early exit(0.2)/water(0.05) / z = torch lightmap
	gl_FragData[4] = vec4(0.0, mix(1.0,0.05,iswater), lmcoord.s, 1.0);
	
}