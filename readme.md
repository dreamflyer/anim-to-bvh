# anim-to-bvh Converter Library

This library provides tools for manipulating `.anim` and `.bvh` files, primarily for **Second Life**. It allows converting animations while preserving Bento bones and ensuring compatibility with **Avastar 3.1.3**, which is currently the most reliable option for importing BVH files into Blender.

## 🚀 Online Tool  
A simple web interface for converting `.anim` files to `.bvh` is available here:  
🔗 **[anim-to-bvh Online Converter](https://dreamflyer.github.io/anim-to-bvh)**  

## 📌 Features  
- **Supports Bento bones** – ensuring compatibility with modern Second Life avatars.  
- **Rotation conversion fully tested** – positional animation *should* work, but has not yet been verified.  
- **Optimized for Avastar 3.1.3** – later versions may have issues with BVH import.  

## ⚠️ Important Notes on Offsets  
`.anim` files do not contain offset information, which can cause issues when importing the resulting `.bvh` into Avastar. (For Second Life itself, offsets are not necessary.)  

To ensure correct offsets in Blender, we recommend using a donor BVH:  
1. Export any animation of your avatar, ensuring that **all bones are included**.  
2. If needed, simply copy a static T-pose and paste it into multiple frames.  
3. Use this donor BVH when converting `.anim` files to apply correct offsets.  

Alternatively, if your avatar is close to the default Second Life skeleton, you can select **"Use default offsets"** in the converter. In this case, no additional setup is needed.  

---

This library is still in active development, and further testing will be done to improve positional animations and ensure broader compatibility.  
Contributions and feedback are welcome! 🚀  
