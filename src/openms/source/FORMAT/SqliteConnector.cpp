// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SqliteConnector.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  void SqliteConnector::openDatabase(const String& filename)
  {
    // Open database
    int rc = sqlite3_open(filename.c_str(), &db_);
    if (rc)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db_));
    }
  }

  bool SqliteConnector::columnExists(sqlite3 *db, const String& tablename, const String& colname)
  {
    bool found = false;

    sqlite3_stmt * xcntstmt;
    String s = "PRAGMA table_info(" + tablename + ")";
    int rc = sqlite3_prepare_v2(db, s.c_str(), -1, &xcntstmt, nullptr);

    if (rc != SQLITE_OK)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Error when determining columns of table: ") + tablename);
    }

    // Go through all columns and check whether the required column exists
    sqlite3_step( xcntstmt );
    while (sqlite3_column_type( xcntstmt, 0 ) != SQLITE_NULL)
    {
      String name = String(reinterpret_cast<const char*>(sqlite3_column_text( xcntstmt, 1 )));
      if (colname == name)
      {
        found = true;
      }

      sqlite3_step( xcntstmt );
    }
    sqlite3_finalize(xcntstmt);
    return found;
  }

  void SqliteConnector::executeStatement(sqlite3 *db, const std::stringstream& statement)
  {
    char *zErrMsg = nullptr;
    std::string insert_str = statement.str();
    int rc = sqlite3_exec(db, insert_str.c_str(), nullptr /* callback */, nullptr, &zErrMsg);
    if (rc != SQLITE_OK)
    {
      String error (zErrMsg);
      std::cerr << "Error message after sqlite3_exec" << std::endl;
      std::cerr << "Prepared statement " << statement.str() << std::endl;
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
    }
  }

  void SqliteConnector::executeStatement(sqlite3 *db, const String& statement)
  {
    char *zErrMsg = nullptr;
    int rc = sqlite3_exec(db, statement.c_str(), nullptr /* callback */, nullptr, &zErrMsg);
    if (rc != SQLITE_OK)
    {
      String error (zErrMsg);
      std::cerr << "Error message after sqlite3_exec" << std::endl;
      std::cerr << "Prepared statement " << statement << std::endl;
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
    }
  }

  void SqliteConnector::executePreparedStatement(sqlite3 *db, sqlite3_stmt** stmt, const String& prepare_statement)
  {
    int rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), stmt, nullptr);
    if (rc != SQLITE_OK)
    {
      std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
      std::cerr << "Prepared statement " << prepare_statement << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
    }
  }

  void SqliteConnector::executeBindStatement(sqlite3 *db, const String& prepare_statement, const std::vector<String>& data)
  {
    sqlite3_stmt *stmt = nullptr;
    const char *curr_loc;
    int rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), &stmt, &curr_loc);
    if (rc != SQLITE_OK)
    {
      std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
      std::cerr << "Prepared statement " << prepare_statement << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
    }

    for (Size k = 0; k < data.size(); k++)
    {
      // Fifth argument is a destructor for the blob.
      // SQLITE_STATIC because the statement is finalized
      // before the buffer is freed:
      rc = sqlite3_bind_blob(stmt, k+1, data[k].c_str(), data[k].size(), SQLITE_STATIC);
      if (rc != SQLITE_OK)
      {
        std::cerr << "SQL error after sqlite3_bind_blob at iteration " << k << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      } 
    }

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE)
    {
      std::cerr << "SQL error after sqlite3_step" << std::endl;
      std::cerr << "Prepared statement " << prepare_statement << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
    }

    // free memory
    sqlite3_finalize(stmt);
  }

}


