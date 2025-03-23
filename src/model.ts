export interface Vector3 {
    x: number;
	y: number;
	z: number;
}

export interface AnimKey {
	x: number,
	y: number,
	z: number,
	time: number
}

export interface AnimJoint {
	joint_name: string;
	joint_priority: number;
	position_keys: AnimKey[],
	rotation_keys: AnimKey[]
}

export interface AnimData {
	constraints: any[];
	duration: number;
	emote_name: string;
	joints: AnimJoint[];
	loop: number;
	sub_version: number;
	version: number;
}

export interface BVHFrame {
	position: Vector3;
	rotation: Vector3;
	time?: number;
}

export interface BVHNode {
	bvhName: string;
	channels?: string[]; 
	bvhTimes?: number[];
	children?: BVHNode[];
	offset: Vector3;
	animFrames?: BVHFrame[];
	bvhFrames: BVHFrame[]
}