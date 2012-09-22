// dbconsole.cpp: главный файл проекта.

#include "stdafx.h"
#include <stdio.h>

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "mysql.h"
/*#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>

#include "mysql_connection.h"

#define EXAMPLE_HOST "localhost"
#define EXAMPLE_USER "root"
#define EXAMPLE_PASS ""
#define EXAMPLE_DB "sprosus"

using namespace System;
using namespace std;
using namespace sql;*/

using namespace System;
using namespace std;

MYSQL mysql; 
MYSQL_RES *res; 
MYSQL_ROW row; 

void exiterr(int exitcode) 
{ 
		fprintf(stderr, "%s\n", mysql_error(&mysql)); 
		exit(exitcode); 
} 

int main()
{
	if (!(mysql_connect(&mysql,"localhost","root",""))) 
			exiterr(1);
    /*string url(argc >= 2 ? argv[1] : EXAMPLE_HOST);
    const string user(argc >= 3 ? argv[2] : EXAMPLE_USER);
    const string pass(argc >= 4 ? argv[3] : EXAMPLE_PASS);
    const string database(argc >= 5 ? argv[4] : EXAMPLE_DB);

	cout << "My first experience with Mysql and C++ " << endl;
    cout << endl;
	
	try {
	
	/* INSERT TUTORIAL CODE HERE! */

		//sql::Driver* driver = get_driver_instance();
		//std::auto_ptr<sql::Connection> con(driver->connect(url, user, pass));
		//con->setSchema(database);
		//std::auto_ptr<sql::Statement> stmt(con->createStatement());

		//stmt->execute("CALL sp_region_get_tree()");
		//std::auto_ptr< sql::ResultSet > res;
		/*do {
				res.reset(stmt->getResultSet());
				while (res->next()) {
					cout << "Result: " << res->getString(1) << endl;
				}
			} while (stmt->getMoreResults());*/
	
   /* }

	catch (sql::SQLException &e) {
     
        cout << "# ERR: SQLException in " << __FILE__;
        cout << "(" << __FUNCTION__ << ") on line " << __LINE__ << endl;
        /* Use what() (derived from std::runtime_error) to fetch the error message */
        //cout << "# ERR: " << e.what();
        //cout << " (MySQL error code: " << e.getErrorCode();
       //s cout << ", SQLState: " << e.getSQLState() << " )" << endl;
	
       /* return EXIT_FAILURE;
    }*/
	
    cout << "Done." << endl;
    
	return true;

}
