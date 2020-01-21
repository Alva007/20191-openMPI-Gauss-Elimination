const fs = require('fs');
const out = process.env.PROG || "inverse";
const child = require('child_process');

const bcd = 120;
const N = [bcd*6, bcd*20, bcd*40];
const P = [2]
// const N = [bcd*6];
// const P = [0];
N.map(n => {
	let works = P.map(p => `mpirun --oversubscribe -np ${p} --hostfile mpi_hosts ./target/cycle-inverse --test ${n}`);
	path = `./test/${out}_${n}.csv`;
	fs.writeFileSync(path, `N: ${n}\np,total,comm\n`, {flag: 'as'});

	works.map(work => {
			let log = child.execSync(work);
			fs.writeFileSync(path, log, {flag: 'as'});
		}
	);

	console.log("DONE: ", n);
})
