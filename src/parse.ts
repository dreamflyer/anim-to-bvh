import {AnimData, BVHNode} from "./model";
import {toEulers, distributeValue} from "./utils";

export function parseAnim(arrayBuffer: ArrayBuffer): AnimData {
    const view = new DataView(arrayBuffer);
    let offset = 0;

    function readU16() { 
        const value = view.getUint16(offset, true);
        offset += 2;
        return value;
    }

    function readU32() { 
        const value = view.getUint16(offset, true);
        offset += 4;
        return value;
    }

    function readS32() { 
        const value = view.getInt32(offset, true);
        offset += 4;
        return value;
    }

    function readF32() { 
        const value = view.getFloat32(offset, true);
        offset += 4;
        return value;
    }

    function readString() {
        let str = "";
        while (offset < view.byteLength) {
            const char = view.getUint8(offset++);
            if (char === 0) break; // NULL-terminated
            str += String.fromCharCode(char);
        }
        return str;
    }

    function readFixedString(size: number) {
        let str = "";
        for (let i = 0; i < size; i++) {
            const char = view.getUint8(offset++);
            if (char !== 0) str += String.fromCharCode(char);
        }
        return str;
    }
	
    const version = readU16();
    const sub_version = readU16();
    const base_priority = readS32();
    const duration = readF32();
    const emote_name = readString();
    const loop_in_point = readF32();
    const loop_out_point = readF32();
    const loop = readS32();
    const ease_in_duration = readF32();
    const ease_out_duration = readF32();
    const hand_pose = readU32();
    const num_joints = readU32();
	
    const joints = [];
	
    for (let i = 0; i < num_joints; i++) {
        const joint_name = readString();
        const joint_priority = readS32();

        const num_rot_keys = readS32();
        const rotation_keys = [];
        for (let j = 0; j < num_rot_keys; j++) {
            const time = readU16();
            const rot_x = readU16();
            const rot_y = readU16();
            const rot_z = readU16();
            rotation_keys.push({time: intToFloat(time, 0, duration), x: intToFloat(rot_x, -1, 1), y: intToFloat(rot_y, -1, 1), z: intToFloat(rot_z, -1, 1)});
        }

        const num_pos_keys = readS32();
        const position_keys = [];
        for (let j = 0; j < num_pos_keys; j++) {
            const time = readU16();
            const pos_x = readU16();
            const pos_y = readU16();
            const pos_z = readU16();
            position_keys.push({time: intToFloat(time, 0, duration), x: intToFloat(pos_x, -5, 5), y: intToFloat(pos_y, -5, 5), z: intToFloat(pos_z, -5, 5)});
        }

        joints.push({ joint_name, joint_priority, rotation_keys, position_keys });
    }
	
    const num_constraints = readS32();
    const constraints = [];
    for (let i = 0; i < num_constraints; i++) {
        const chain_length = view.getUint8(offset++);
        const constraint_type = view.getUint8(offset++);
        const source_volume = readFixedString(16);
        const source_offset = [readF32(), readF32(), readF32()];
        const target_volume = readFixedString(16);
        const target_offset = [readF32(), readF32(), readF32()];
        const target_dir = [readF32(), readF32(), readF32()];
        const ease_in_start = readF32();
        const ease_in_stop = readF32();
        const ease_out_start = readF32();
        const ease_out_stop = readF32();

        constraints.push({
            chain_length, constraint_type, source_volume, source_offset,
            target_volume, target_offset, target_dir, ease_in_start, ease_in_stop,
            ease_out_start, ease_out_stop
        });
    }
	
	joints.forEach((item: any) => item.rotation_keys.forEach((rot: any) => {
		if(!item.euler_keys) {
			item.euler_keys = [];
		}
		
		item.euler_keys.push(toEulers(rot));
	}));

    return { version, sub_version, duration, emote_name, loop, joints, constraints };
}


function intToFloat(val: number, min: number, max: number): number {
	const one = (max - min) / 65535.0;
	
	const result =  min + val * one;
	
	if(Math.abs(result) < one) {
		return 0;
	}
	
	return result;
}

function enumerate(content: string, key: string, alter: string): string {
	let result: string = content;
	
	let count = 0;
	
	while(result.includes("\"" + key + "\"")) {
			result = result.replace("\"" + key + "\"", "\"" + alter + count + "\"");
			
			count += 1;
	}
	
	return result;
}

function compose(node: any, name: string): any {
	const children: any[] = [];
	const result: any = {};
	
	let cnts: string[] = [];
	let jnts: string[] = [];
	
	if(Object.keys(node).includes("End Site")) {
		node.jnt1 = "end";
	}
	
	Object.keys(node).forEach(item => {		
		if(item.includes("cnt")) {
			cnts.push(item);
			
			return;
		}
		
		if(item.includes("jnt")) {
			jnts.push(item);
		}
	});
	
	cnts = cnts.sort((val1, val2) => {return parseInt(val1.replace("cnt", "")) - parseInt(val2.replace("cnt", ""));});
	jnts = jnts.sort((val1, val2) => {return parseInt(val1.replace("jnt", "")) - parseInt(val2.replace("jnt", ""));});
	
	jnts.forEach((item, i) => children.push(compose(node[cnts[i]], node[item].trim())));
	
	if(node.OFFSET) {
		const offset: any = node.OFFSET.trim().split(" ").filter((item: any) => !isNaN(parseFloat(item))).map(parseFloat);
		
		result.offset = {x: offset[0], y: offset[1], z: offset[2]};
	}
	
	if(node.CHANNELS) {
		const channels: any = node.CHANNELS.trim().split(" ").map((item: any) => item.trim()).filter((item: any) => !!item).filter((item: any) => isNaN(parseFloat(item)));
		
		result.channels = channels;
	}
	
	result.bvhName = name;
	
	if(children.length) {
		result.children = children;
	}
	
	return result;
}

function cleanup(data: any) {	
	data.jnt0 = data.ROOT;
	
	return compose(data, "root");
}

function parseFrames(rows: string[]): number[][] {
	const splitedRows: string[][] = rows.map(item => item.split(" ").map(item => item.trim()).filter(item => !!item));
	
	return splitedRows.map(item => item.map(parseFloat));
}

function parseFramesPart(framesPart: string): {framesLength: number, frameDuration: number, frames: number[][]} {
	const framesRows: string[] = framesPart.split("\n");
	
	let timeIndex: number = -1;
	
	for(let i = 0; i < framesRows.length; i++) {
		if(framesRows[i].toLowerCase().includes("time")) {
			timeIndex = i;
			
			break;
		}
	}
	
	if(timeIndex < 0) {
		return {framesLength: 0, frameDuration: 0, frames: []};
	}
	
	const framesLength: number = <number>parseInt(framesRows[timeIndex - 1].split(" ").map((item: string) => item.trim()).filter((item: string) => !!item).filter((item: string) => !isNaN(<number><any>item))[0]);
	const frameDuration: number = <number>parseFloat(framesRows[timeIndex].split(" ").map((item: string) => item.trim()).filter((item: string) => !!item).filter((item: string) => !isNaN(<number><any>item))[0]);
	
	while(!framesRows[0].toLowerCase().includes("time")) {
		framesRows.shift();
	}
	
	framesRows.shift();
	
	const frames: number[][] = parseFrames(framesRows);
	
	return {framesLength, frameDuration, frames}
}

export function distributeSingleFrame(hierarchy: BVHNode, frame: number[]) {
	hierarchy.children?.toReversed().forEach((child: any) => distributeSingleFrame(child, frame));
	
	const position = {x: 0, y: 0, z: 0};
	const rotation = {x: 0, y: 0, z: 0};
	
	hierarchy.channels?.toReversed().forEach((item: any) => {
		const value: number = <number><any>frame.pop();
		
		distributeValue(position, rotation, item, value);
	});
	
	if(!hierarchy.bvhFrames) {
		hierarchy.bvhFrames = [];
	}
	
	hierarchy.bvhFrames.push({position, rotation});
}

export function parseBVH(text: string): BVHNode {
	const parts: string[] = text.replaceAll("\r", "").replaceAll("\t", " ").split("MOTION");
	
	let result = parts[0].split("HIERARCHY")[1];
	
	"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 ".split("").forEach(item => {
		result = result.replaceAll(item + "\n", item + "\",\n");
	});
	
	result = result.replaceAll("JOINT","\"JOINT\":\"");
	result = result.replaceAll("OFFSET","\"OFFSET\":\"");
	result = result.replaceAll("CHANNELS","\"CHANNELS\":\"");
	result = result.replaceAll("ROOT","\"ROOT\":\"");
	result = result.replaceAll("End Site","\"End Site\":\"");
	result = result.replaceAll("{","\"content\": {");
		
	result = result.split("}").map(item => {
		if(item.trim().endsWith(",")) {
			return item.trim() + "\"dummy\": {}";
		}
		
		return item;
	}).join("}");
	
	result = result.split("}").map(item => {
		if(item.trim().startsWith("\"JOINT\"")) {
			return "," + item.trim();
		}
		
		return item;
	}).join("}");
	
	let count = 0;
	
	result = enumerate(result, "JOINT", "jnt");
	result = enumerate(result, "content", "cnt");
	
	const hierarchy: BVHNode = cleanup(JSON.parse("{" + result + "}")).children[0];
	
	const animation = parseFramesPart(parts[1]);
	
	hierarchy.bvhTimes = [];
	
	for(let i = 0; i < animation.framesLength; i++) {
		hierarchy.bvhTimes.push(animation.frameDuration * i);
	}
	
	animation.frames.forEach(item => distributeSingleFrame(hierarchy, item));
		
	return hierarchy;
}