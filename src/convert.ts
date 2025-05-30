import {AnimData, AnimJoint, AnimKey, BVHNode, BVHFrame, Vector3} from "./model";
import {append, toQuaternion, lerpValues, getUniformTimes, clipTimesToClosestBVHTime, lerpVector, lerpQuaternion, quaternionToEulers, floatToString} from "./utils";

import {hierarchy} from "./hierarchy";
import {aliases} from "./aliases";

import {Quaternion} from "quaternion";

function offsetToString(offset: Vector3, digits: number): string {
	return "OFFSET " + floatToString(offset.x, digits) + " " + floatToString(offset.y, digits) + " " + floatToString(offset.z, digits);
}

function channelsString(node: BVHNode): string {
	if(!node.channels) {
		return "";
	}
	
	if(node.bvhName === "hip") {
		return "CHANNELS 6 Xposition Yposition Zposition Xrotation Yrotation Zrotation";
	}
	
	return "CHANNELS " + node.channels.length + " " + node.channels.join(" ");
}

function appendNode(joint: BVHNode, tabs: string): string {
	let result = "";
	
	const boneType = (joint.bvhName === "hip") ? "ROOT" : "JOINT";
	
	const channels = channelsString(joint);
	
	const offset = (joint.bvhName === "hip") ? offsetToString(joint.offset, 6) : offsetToString(joint.offset, 4);
	
	if(joint.bvhName != "end") {
		result += tabs + boneType + " " + joint.bvhName + "\n" + tabs + "{\n";
	} else {
		result += tabs + "End Site" + "\n" + tabs + "{\n";
	}
	
	result += tabs + "\t" + offset + "\n";
	
	if(joint.bvhName != "end") {
		result += tabs + "\t" + channels + "\n";
	}
	
	if(joint.children) {
		joint.children.forEach((item: any) => {result+=appendNode(item, tabs + "\t")});
	}
	
	result += tabs + "}\n";
	
	return result;
}

function containsNames(node: any, bvhNames: string[]) {
	if(bvhNames.includes(node.bvhName)) {
		return true;
	}
	
	if(!node.children) {
		return false;
	}
	
	return !!node.children.map((item: any) => containsNames(item, bvhNames)).find((item: any) => !!item);
}

function collectNodes(node: any, bvhNames: string[]): any {
	const result: any = {};
	
	if(containsNames(node, bvhNames)) {
		result.bvhName = node.bvhName;
	} else {
		result.exclude = true;
		
		return result;
	}
	
	if(node.children && !!node.children.map((item: any) => containsNames(item, bvhNames)).find((item: any) => !!item)) {
		result.children = node.children.map((item: any) => collectNodes(item, bvhNames)).filter((item: any) => !item.exclude);
	} else {
		result.children = [];
	}
	
	if(result.children.length > 0) {
		return result;
	}
	
	result.children.push({bvhName: "end"});
	
	return result;
}

function subTree(joints: any[]): any {
	const names: string[] = joints.map(item => item.joint_name);
	const bvhNames: string[] = names.map(item => aliases[item] || item);
	
	return collectNodes(hierarchy, bvhNames);
}

export function visitNode(node: BVHNode, visitor: (node: BVHNode) => void, childrenFirst: boolean = false): void {
	if(node.children && childrenFirst) {
		node.children.toReversed().forEach((item: any) => visitNode(item, visitor, true));
	}
		
	visitor(node);
	
	if(node.children && !childrenFirst) {
		node.children.toReversed().forEach((item: any) => visitNode(item, visitor, false));
	}
}

function extractFramesLength(animJoints: any): number {
	const joint: any = animJoints.find((item: any) => item.position_keys?.length || item.rotation_keys?.length);
		
	return joint?.position_keys?.length || joint?.rotation_keys?.length;
}

function extractTimes(animJoints: any): number[] {
	const joint: any = animJoints.find((item: any) => item.position_keys?.length || item.rotation_keys?.length);
	
	const timeHolders = joint?.position_keys?.length ? joint.position_keys : joint?.rotation_keys;
	
	return (timeHolders || []).map((item: any) => item.time);
}

function fillChannels(node: any, joint: any): void {
	if(node.bvhName == "hip") {
		node.channels = ["Xposition", "Yposition", "Zposition", "Xrotation", "Yrotation", "Zrotation"];
		
		return;
	}
	
	node.channels = [];
	
	if(joint?.position_keys?.length) {
		node.channels.push("Xposition");
		node.channels.push("Yposition");
		node.channels.push("Zposition");
	}
	
	node.channels.push("Xrotation");
	node.channels.push("Yrotation");
	node.channels.push("Zrotation");
}

function animPositionToBvh(position: Vector3): Vector3 {
	const multiplier = 39.3795;
	
	return {x: position.y * multiplier, y: position.z * multiplier, z: position.x * multiplier}
}

function fillKeyFrames(data: AnimData, bvhNode: BVHNode, fps: number): void {
	const animJoints: any[] = data.joints;
	
	const length: number = extractFramesLength(animJoints);
	
	const bvhTimes: number[] = getUniformTimes(data.duration || 2, 1 / fps);
	
	visitNode(bvhNode, (node) => {
		const joint: AnimJoint = animJoints.find((item: any) => aliases[item.joint_name] === node.bvhName);
		
		node.offset = {x: 0, y: 0, z:0};
		
		if(node.bvhName != "end") {
			fillChannels(node, joint);
			
			node.animKeys = {
				positions: joint?.position_keys || [],
				rotations: joint?.rotation_keys || []
			}
			
			const positionTimes: number[] = clipTimesToClosestBVHTime(node.animKeys.positions.map((item: AnimKey) => item.time), bvhTimes);
			const rotationTimes: number[] = clipTimesToClosestBVHTime(node.animKeys.rotations.map((item: AnimKey) => item.time), bvhTimes);
			
			const positions: Vector3[] = lerpValues(node.animKeys.positions.map((item: any) => animPositionToBvh(item)), positionTimes, bvhTimes, {x: 0, y: 0, z: 0},lerpVector);
			const rotations: Vector3[] = lerpValues(node.animKeys.rotations.map((item: any) => toQuaternion(item)), rotationTimes, bvhTimes, new Quaternion(), lerpQuaternion).map(item => quaternionToEulers(item));
						
			node.bvhFrames = [];
						
			bvhTimes.forEach((item: number, i: number) => node.bvhFrames.push({
				position: positions[i],
				rotation: rotations[i]
			}));
		}
	});
	
	bvhNode.bvhTimes = bvhTimes;
}

function getValue(bvhNode: BVHNode, channel: string, frameNum: number): number {
	const frame = bvhNode.bvhFrames[frameNum];
	
	const key = channel.toLowerCase()[0];
	
	const data = (channel.includes("pos")) ? frame.position : frame.rotation;
	
	const value: number = <number>(<any>data)[key];
	
	return (Math.abs(value) > 0.00000001) ? value : 0;	
}

function getValues(bvhNode: BVHNode, frameNum: number) {
	return bvhNode.channels!.map((item: any) => getValue(bvhNode, item, frameNum));
}

function getFrameValues(bvhNode: any, frameNum: number): number[] {
	const result: number[] = [];
	
	visitNode(bvhNode, (node: any) => {
		if(!node.channels) {
			return;
		}
		
		result.unshift(...getValues(node, frameNum));
	}, true);
	
	return result;
}

function getFrameRow(bvhNode: any, frameNum: number): string {
	const values: number[] = getFrameValues(bvhNode, frameNum);
	
	return values.map(item => floatToString(item, 4)).join(" ") + " \n";
}



export function serializeBVH(bvhNode: BVHNode): string {
	let result = "HIERARCHY\n";
	
	result += appendNode(bvhNode, "");
	result += "MOTION\n";
	result += "Frames: " + bvhNode.bvhTimes!.length + "\n";
	result += "Frame Time: " + floatToString(bvhNode.bvhTimes![1], 6) + "\n";
	
	for(let i = 0; i < bvhNode.bvhTimes!.length; i++) {
		result += getFrameRow(bvhNode, i);
	}
	
	return result;
}

export function toBVH(data: AnimData, fps: number = 24): BVHNode {
	const bvhNode: BVHNode = <BVHNode><any>subTree(data.joints);
		
	fillKeyFrames(data, bvhNode, fps);
		
	return bvhNode;
}

export function collectOffsets(bvhNode: BVHNode): {[name: string]: Vector3} {
	const result: {[name: string]: Vector3} = {};
	
	visitNode(bvhNode, node => {
		if(node.children && node.children.length && (node.children[0].bvhName == "end")) {
			node.children[0].parentName = node.bvhName;
		}
		
		if(node.bvhName == "end") {
			result["end_" + node.parentName] = node.offset;
			
			return;
		}
		
		result[node.bvhName] = node.offset;
	});
	
	return result;
}

export function collectReferenceFrame(bvhNode: BVHNode): {[name: string]: BVHFrame} {
	const result: {[name: string]: BVHFrame} = {};
	
	visitNode(bvhNode, node => {
		if(node.bvhName == "end") {
			return;
		}
		
		const frame: BVHFrame = {
			position: node.bvhFrames?.[0]?.position || {x: 0, y: 0, z: 0},
			rotation: node.bvhFrames?.[0]?.rotation || {x: 0, y: 0, z: 0}
		}
		
		result[node.bvhName] = frame;
	});
	
	return result;
}
