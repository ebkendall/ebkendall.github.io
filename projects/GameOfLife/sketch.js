function make2DArray(cols, rows) {
	let arr = new Array(cols);
	for (let i = 0; i < arr.length; i++) {
		arr[i] = new Array(rows);
	}
	return(arr)
}

let grid;
let next;
let cols;
let rows;
let resolution = 20;
let fillCol = [50, 255];

// Switch
let turnOff;
let rand_or_blank;

// Button selection
let button; 		 // "turnOff" switch
let button_template; // "rand_or_blank" switch



function setup() {
	createCanvas(800, 600);
	cols = width / resolution;
	rows = height / resolution;

	grid = make2DArray(cols, rows);

	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j] = new Cell(i * resolution, j * resolution);
			grid[i][j].val = floor(random(2));
		}
	}

	next = grid;

	turnOff = 1;
	rand_or_blank = 1;

	// button characteristics
	button = createButton('START/STOP');
	button.mousePressed(changeB);

	// button_template = createButton('RANDOM or CLEAR');
	// button_template.mousePressed(changeB2);

}

function draw() {
	background(0);

	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].display();
		}
	}

	if(turnOff == 0) {
		next = make2DArray(cols, rows);

		for (let i = 0; i < cols; i++) {
			for (let j = 0; j < rows; j++) {

				let state = grid[i][j].val;

				let sum = 0;
				let neighbors = countNeighbors(grid, i , j);

				if (state == 0 && neighbors == 3) {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = 1;
				} else if (state == 1 && (neighbors < 2 || neighbors > 3)) {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = 0;
				} else {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = state;
				}
			}
		}

		grid = next;
	}
}

function countNeighbors(grid, x, y) {

	let sum = 0;

	for (let i = -1; i < 2; i ++) {
		for (let j = -1; j < 2; j++) {

			let w = (x + cols + i) % cols;
			let z = (y + rows + j) % rows;

			sum += grid[w][z].val;
		}
	}
	sum -= grid[x][y].val;
	return sum;
}

function changeB() {
	if (turnOff == 0) {
		turnOff = 1;
	} else {
		turnOff = 0;
	}
}

function changeB2() {
	if (rand_or_blank == 0) {
		rand_or_blank = 1;
	} else {
		rand_or_blank = 0;
	}
}

function mousePressed() {
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].clicked(mouseX, mouseY);
		}
	}
}

class Cell {
	constructor(x, y) {
		this.x = x;
		this.y = y;
		this.val = 0;
	}

	clicked(x, y) {
		if((x > this.x) && (x < this.x + resolution) &&
		   (y > this.y) && (y < this.y + resolution)) {
			   if(this.val == 1) {
				   this.val = 0;
			   } else {
				   this.val = 1;
			   }
		}

	}

	display() {
		fill(fillCol[this.val]);
		stroke(0);
		rect(this.x, this.y, resolution-1, resolution-1);
	}
}

// I want to add the following:
// (1) Menu for starting screens (only appears during STOP)
// (2) Blank template
// (3) Random template
