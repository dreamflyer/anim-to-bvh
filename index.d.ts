// Generated by dts-bundle-generator v9.5.1

export interface Vector3 {
	x: number;
	y: number;
	z: number;
}
export interface AnimKey {
	x: number;
	y: number;
	z: number;
	time: number;
}
export interface AnimJoint {
	joint_name: string;
	joint_priority: number;
	position_keys: AnimKey[];
	rotation_keys: AnimKey[];
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
export interface AnimKeys {
	positions: AnimKey[];
	rotations: AnimKey[];
}
export interface BVHNode {
	bvhName: string;
	channels?: string[];
	bvhTimes?: number[];
	children?: BVHNode[];
	offset: Vector3;
	animKeys?: AnimKeys;
	bvhFrames: BVHFrame[];
}
export declare function parseAnim(arrayBuffer: ArrayBuffer): AnimData;
export declare function parseBVH(text: string): BVHNode;
export declare function visitNode(node: BVHNode, visitor: (node: BVHNode) => void, childrenFirst?: boolean): void;
export declare function serializeBVH(bvhNode: BVHNode): string;
export declare function toBVH(data: AnimData, fps?: number): BVHNode;
export declare const femaleOffsets: {
	[key: string]: Vector3;
};
export declare const maleOffsets: {
	[key: string]: Vector3;
};

export {};
