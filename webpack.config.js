const path = require("path");

module.exports = {
  entry: "./src/index.ts",

  output: {
    filename: "bundle.js",
    path: path.resolve(__dirname, "dist"),
    publicPath: "/",
	library: "AnimToBvh",
	globalObject: "typeof self !== 'undefined' ? self : this",
	libraryTarget: "umd"
  },

  mode: "development",
  
  resolve: {
    extensions: [".ts", ".js"],
  },
  
  module: {
    rules: [
      {
        test: /\.ts$/,
        use: "ts-loader",
        exclude: /node_modules/,
      },
    ],
  },

  devtool: "inline-source-map",
};